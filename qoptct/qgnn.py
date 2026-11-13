"""Сеть, не зависящая от числа кубитов (графовая нейросеть).

Вход — то же наблюдение, что у обычной сети: вещественная и мнимая части 2^n амплитуд.
Внутри сеть сама считает «локальные» свойства состояния:
  - для каждого кубита — вектор Блоха (средние X, Y, Z), 3 числа;
  - для каждой пары кубитов — корреляции <A_i B_j> для A,B из {X,Y,Z}, 9 чисел.
Узлы графа = кубиты, рёбра = пары. Несколько шагов «обмена сообщениями» с ОБЩИМИ весами,
затем:
  - ход «одиночная команда на кубите q» оценивается по признакам кубита q;
  - ход «CX(c, t)» — по признакам c, t и ребра c->t.
Так один и тот же набор весов работает для 5, 6, 7, 8, ... кубитов, а перестановка кубитов
меняет ответ ровно так же, как меняются сами ходы (симметрия встроена).
Порядок выходов совпадает с qfast.Engine.gates.
"""

import itertools

import torch
import torch.nn as nn

PAULI = torch.tensor([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], dtype=torch.complex64)
PP = torch.stack([torch.kron(a, b) for a in PAULI for b in PAULI])  # (9, 4, 4)


def mlp(i, h, o):
    return nn.Sequential(nn.Linear(i, h), nn.ReLU(), nn.Linear(h, h), nn.ReLU(), nn.Linear(h, o))


def local_features(x, n):
    """x: (B, 2*2^n) -> node (B, n, 4), edge (B, n, n, 9) (edge[c, t] — корреляции в порядке (c, t))."""
    D = 2 ** n
    z = torch.complex(x[:, :D], x[:, D:])
    B = z.shape[0]
    psi = z.reshape([B] + [2] * n)  # ось 1+(n-1-q) — кубит q (порядок Qiskit)
    ax = lambda q: 1 + (n - 1 - q)
    dev = x.device
    P1, P2 = PAULI.to(dev), PP.to(dev)
    node = torch.zeros(B, n, 4, device=dev)
    for q in range(n):
        m = psi.movedim(ax(q), 1).reshape(B, 2, -1)
        rho = m @ m.conj().transpose(1, 2)
        bloch = torch.einsum("bij,kji->bk", rho, P1).real
        node[:, q, :3] = bloch
        node[:, q, 3] = (rho @ rho).diagonal(dim1=1, dim2=2).sum(-1).real  # чистота кубита
    edge = torch.zeros(B, n, n, 9, device=dev)
    for i, j in itertools.combinations(range(n), 2):
        m = psi.movedim((ax(i), ax(j)), (1, 2)).reshape(B, 4, -1)
        rho = m @ m.conj().transpose(1, 2)
        c = torch.einsum("bij,kji->bk", rho, P2).real  # <A_i B_j>, A — первый индекс
        edge[:, i, j] = c
        edge[:, j, i] = c.view(B, 3, 3).transpose(1, 2).reshape(B, 9)
    return node, edge


class QGNN(nn.Module):
    arch = "gnn"

    def __init__(self, d=128, layers=3):
        super().__init__()
        self.d, self.layers_n = d, layers
        self.node_in = mlp(4, d, d)
        self.edge_in = mlp(9, d, d)
        self.msg = nn.ModuleList([mlp(3 * d, d, d) for _ in range(layers)])
        self.upd = nn.ModuleList([mlp(2 * d, d, d) for _ in range(layers)])
        self.head1 = mlp(2 * d, d, 6)      # h, x, s, sdg, t, tdg на кубите q
        self.headcx = mlp(4 * d, d, 1)     # CX(c, t)

    def forward(self, x):
        D2 = x.shape[1]
        n = (D2 // 2).bit_length() - 1
        node, edge = local_features(x, n)
        B = x.shape[0]
        h = self.node_in(node)                      # (B, n, d)
        e = self.edge_in(edge)                      # (B, n, n, d)
        mask = (1 - torch.eye(n, device=x.device)).view(1, n, n, 1)
        for L in range(self.layers_n):
            hi = h.unsqueeze(2).expand(B, n, n, self.d)
            hj = h.unsqueeze(1).expand(B, n, n, self.d)
            m = self.msg[L](torch.cat([hi, hj, e], -1)) * mask
            agg = m.sum(2) / max(1, n - 1)
            h = h + self.upd[L](torch.cat([h, agg], -1))
        g = h.mean(1, keepdim=True).expand(B, n, self.d)
        one = self.head1(torch.cat([h, g], -1)).reshape(B, n * 6)   # порядок: кубит, затем команда
        cs, ts = zip(*[(c, t) for c in range(n) for t in range(n) if c != t])
        cs, ts = list(cs), list(ts)
        gc = h.mean(1)
        cx = self.headcx(torch.cat([h[:, cs], h[:, ts], e[:, cs, ts], gc.unsqueeze(1).expand(B, len(cs), self.d)],
                                   -1)).squeeze(-1)
        return torch.cat([one, cx], 1)


def selftest():
    """Проверка порядка выходов и симметрии перестановок на случайных состояниях."""
    import numpy as np
    from qfast import Engine
    from ei2 import obs_batch, perm_maps
    torch.manual_seed(0)
    net = QGNN(d=32, layers=2).eval()
    for n in (3, 5, 7):
        eng = Engine(n)
        rng = np.random.default_rng(n)
        s = rng.normal(size=(4, 2 ** n)) + 1j * rng.normal(size=(4, 2 ** n))
        s /= np.linalg.norm(s, axis=1, keepdims=True)
        perm = rng.permutation(n)
        jmap, amap = perm_maps(eng, perm)
        sp = np.empty_like(s); sp[:, jmap] = s
        with torch.no_grad():
            a = net(torch.from_numpy(obs_batch(s)))
            b = net(torch.from_numpy(obs_batch(sp)))
        assert a.shape[1] == eng.n_actions
        err = (b[:, torch.from_numpy(amap)] - a).abs().max().item()
        print(f"{n} кубитов: выходов {a.shape[1]} = действий {eng.n_actions}; ошибка симметрии {err:.1e}")


if __name__ == "__main__":
    selftest()


class QGNN2(nn.Module):
    """Облегчённая версия (этап 6б): тот же смысл, но обмен сообщениями дешевле.
    - признаки рёбер переводятся в вектора ОДИН раз для всех слоёв (одна матричная операция);
    - сообщение = ReLU(A·h_i + B·h_j + C_l·e_ij), затем одна линейная операция — вместо трёхслойной
      сети на каждую пару кубитов на каждом шаге;
    - узкие рёбра (de < d).
    """
    arch = "gnn2"
    features = staticmethod(local_features)

    def __init__(self, d=128, de=64, layers=3):
        super().__init__()
        self.d, self.de, self.layers_n = d, de, layers
        self.node_in = mlp(4, d, d)
        self.edge_in = nn.Sequential(nn.Linear(9, de), nn.ReLU(), nn.Linear(de, de))
        self.edge_proj = nn.Linear(de, layers * d)      # C_l e для всех слоёв сразу
        self.a = nn.ModuleList([nn.Linear(d, d, bias=False) for _ in range(layers)])
        self.b = nn.ModuleList([nn.Linear(d, d) for _ in range(layers)])
        self.out = nn.ModuleList([nn.Linear(d, d) for _ in range(layers)])
        self.upd = nn.ModuleList([mlp(2 * d, d, d) for _ in range(layers)])
        self.head1 = mlp(2 * d, d, 6)
        self.cx_c = nn.Linear(d, d, bias=False)
        self.cx_t = nn.Linear(d, d, bias=False)
        self.cx_e = nn.Linear(de, d)
        self.cx_g = nn.Linear(d, d, bias=False)
        self.cx_out = nn.Sequential(nn.ReLU(), nn.Linear(d, d), nn.ReLU(), nn.Linear(d, 1))

    def forward(self, x):
        n = (x.shape[1] // 2).bit_length() - 1
        node, edge = self.features(x, n)
        B, d = x.shape[0], self.d
        h = self.node_in(node)
        e = self.edge_in(edge)                                    # (B, n, n, de)
        ce = self.edge_proj(e).view(B, n, n, self.layers_n, d)    # один раз
        mask = (1 - torch.eye(n, device=x.device)).view(1, n, n, 1)
        for L in range(self.layers_n):
            m = torch.relu(self.a[L](h).unsqueeze(2) + self.b[L](h).unsqueeze(1) + ce[:, :, :, L]) * mask
            agg = self.out[L](m.sum(2) / max(1, n - 1))
            h = h + self.upd[L](torch.cat([h, agg], -1))
        g = h.mean(1)
        one = self.head1(torch.cat([h, g.unsqueeze(1).expand(B, n, d)], -1)).reshape(B, n * 6)
        cs, ts = zip(*[(c, t) for c in range(n) for t in range(n) if c != t])
        cs, ts = list(cs), list(ts)
        z = self.cx_c(h)[:, cs] + self.cx_t(h)[:, ts] + self.cx_e(e[:, cs, ts]) + self.cx_g(g).unsqueeze(1)
        return torch.cat([one, self.cx_out(z).squeeze(-1)], 1)


_IDX_CACHE = {}


def local_features_fast(x, n):
    """То же, что local_features, но все пары кубитов — матричными операциями в ВЕЩЕСТВЕННЫХ числах.
    Для каждого кубита q: X_q·ψ (переворот оси q), Y_q·ψ, Z_q·ψ (знаки). Тогда
    <A_i B_j> = Re<A_i ψ | B_j ψ> = Re·Re + Im·Im — одна пара матричных умножений вместо цикла по парам.
    Вектор Блоха: <A_q> = Re<A_q ψ | ψ>; чистота кубита = (1 + |r|^2) / 2. (Этап 6б/10в; в 1,65 раза быстрее
    прежней комплексной версии, совпадение до 1e-7.)"""
    D = 2 ** n
    B = x.shape[0]
    dev = x.device
    key = (n, str(dev))
    if key not in _IDX_CACHE:
        k = torch.arange(D, device=dev)
        _IDX_CACHE[key] = torch.stack([1 - 2 * ((k >> q) & 1) for q in range(n)]).float()
    sign = _IDX_CACHE[key]
    re, im = x[:, :D], x[:, D:]
    pre, pim = re.reshape([B] + [2] * n), im.reshape([B] + [2] * n)
    xr = torch.stack([pre.flip(1 + (n - 1 - q)).reshape(B, D) for q in range(n)], 1)
    xi = torch.stack([pim.flip(1 + (n - 1 - q)).reshape(B, D) for q in range(n)], 1)
    zr, zi = re.unsqueeze(1) * sign, im.unsqueeze(1) * sign
    yr, yi = xi * sign, -xr * sign
    Pr = torch.stack([xr, yr, zr], 2).reshape(B, 3 * n, D)
    Pi = torch.stack([xi, yi, zi], 2).reshape(B, 3 * n, D)
    G = Pr @ Pr.transpose(1, 2) + Pi @ Pi.transpose(1, 2)
    bloch = (Pr @ re.unsqueeze(2) + Pi @ im.unsqueeze(2)).view(B, n, 3)
    node = torch.cat([bloch, ((1 + (bloch ** 2).sum(-1)) / 2).unsqueeze(-1)], -1)
    edge = G.view(B, n, 3, n, 3).permute(0, 1, 3, 2, 4).reshape(B, n, n, 9)
    edge = edge * (1 - torch.eye(n, device=dev)).view(1, n, n, 1)
    return node, edge


class QGNN3(QGNN2):
    """QGNN2 + быстрый расчёт свойств (local_features_fast). Совместима по весам с QGNN2."""
    arch = "gnn3"

    @staticmethod
    def features(x, n):
        # на 9+ кубитах быстрый способ в ~2 раза быстрее; на малых — обычный
        return local_features_fast(x, n) if n >= 9 else local_features(x, n)


class QGNNU(QGNN3):
    """Для синтеза ОПЕРАЦИЙ (этап 10): вход — состояние 2n кубитов (n рабочих + n двойников).
    Узлам добавлен признак «рабочий / двойник»; выходы — только команды на рабочих кубитах,
    в порядке Engine(n).gates. Веса общие → одна сеть на любое n.
    (Этап 10в: выходы считаются только для рабочих кубитов и их пар — без лишних вычислений.)"""
    arch = "gnnu"

    def __init__(self, d=128, de=64, layers=3):
        super().__init__(d=d, de=de, layers=layers)
        self.node_in = mlp(5, d, d)

    @staticmethod
    def features(x, N):
        node, edge = local_features_fast(x, N) if N >= 7 else local_features(x, N)
        flag = (torch.arange(N, device=x.device) < N // 2).float().view(1, N, 1).expand(node.shape[0], N, 1)
        return torch.cat([node, flag], -1), edge

    def forward(self, x):
        N = (x.shape[1] // 2).bit_length() - 1
        n = N // 2
        node, edge = self.features(x, N)
        B, d = x.shape[0], self.d
        h = self.node_in(node)
        e = self.edge_in(edge)
        ce = self.edge_proj(e).view(B, N, N, self.layers_n, d)
        mask = (1 - torch.eye(N, device=x.device)).view(1, N, N, 1)
        for L in range(self.layers_n):
            m = torch.relu(self.a[L](h).unsqueeze(2) + self.b[L](h).unsqueeze(1) + ce[:, :, :, L]) * mask
            agg = self.out[L](m.sum(2) / max(1, N - 1))
            h = h + self.upd[L](torch.cat([h, agg], -1))
        g = h.mean(1)
        hs = h[:, :n]
        one = self.head1(torch.cat([hs, g.unsqueeze(1).expand(B, n, d)], -1)).reshape(B, n * 6)
        cs, ts = zip(*[(c, t) for c in range(n) for t in range(n) if c != t])
        cs, ts = list(cs), list(ts)
        z = self.cx_c(h)[:, cs] + self.cx_t(h)[:, ts] + self.cx_e(e[:, cs, ts]) + self.cx_g(g).unsqueeze(1)
        return torch.cat([one, self.cx_out(z).squeeze(-1)], 1)


class QGNNH(QGNNU):
    arch="gnnh"
    def __init__(self,d=128,de=64,layers=3,gw=512,ns=(3,4,5)):
        super().__init__(d=d,de=de,layers=layers)
        self.glob=nn.ModuleDict({str(n):nn.Sequential(nn.Linear(2*4**n,gw),nn.ReLU(),nn.Linear(gw,gw),nn.ReLU()) for n in ns})
        self.g2node=nn.Linear(gw,2*d); nn.init.zeros_(self.g2node.weight); nn.init.zeros_(self.g2node.bias)
    def forward(self,x):
        N=(x.shape[1]//2).bit_length()-1; n=N//2
        node,edge=self.features(x,N); B,d=x.shape[0],self.d
        gv=self.g2node(self.glob[str(n)](x)).view(B,1,2,d)
        flag=(torch.arange(N)<n).view(1,N,1)
        h=self.node_in(node)+torch.where(flag,gv[:,:,0],gv[:,:,1])+(self._pos if hasattr(self,'_pos') else 0)   # рабочие / двойники получают свою часть
        e=self.edge_in(edge); ce=self.edge_proj(e).view(B,N,N,self.layers_n,d)
        mask=(1-torch.eye(N)).view(1,N,N,1)
        for L in range(self.layers_n):
            m=torch.relu(self.a[L](h).unsqueeze(2)+self.b[L](h).unsqueeze(1)+ce[:,:,:,L])*mask
            h=h+self.upd[L](torch.cat([h,self.out[L](m.sum(2)/max(1,N-1))],-1))
        g=h.mean(1); hs=h[:,:n]
        one=self.head1(torch.cat([hs,g.unsqueeze(1).expand(B,n,d)],-1)).reshape(B,n*6)
        cs,ts=zip(*[(c,t) for c in range(n) for t in range(n) if c!=t]); cs,ts=list(cs),list(ts)
        z=self.cx_c(h)[:,cs]+self.cx_t(h)[:,ts]+self.cx_e(e[:,cs,ts])+self.cx_g(g).unsqueeze(1)
        return torch.cat([one,self.cx_out(z).squeeze(-1)],1)

class QGNNHP(QGNNH):
    """Глобальная ветвь выдаёт СВОЙ вектор каждому кубиту (по позиции) — не эквивариантна; симметрию даём
    аугментацией перестановок (ei2 --perm-aug). Нулевая инициализация — старт = QGNNU."""
    arch="gnnhp"
    def __init__(self,d=128,de=64,layers=3,gw=512,ns=(3,4,5)):
        super().__init__(d=d,de=de,layers=layers,gw=gw,ns=ns)
        self.g2pos=nn.ModuleDict({str(n):nn.Linear(gw,2*n*d) for n in ns})
        for m in self.g2pos.values(): nn.init.zeros_(m.weight); nn.init.zeros_(m.bias)
    def forward(self,x):
        N=(x.shape[1]//2).bit_length()-1; n=N//2; B,d=x.shape[0],self.d
        self._pos=self.g2pos[str(n)](self.glob[str(n)](x)).view(B,N,d)
        return super().forward(x)


class QValueU(QGNNU):
    """Оценщик перспективности (этап 10г): то же «тело», что у QGNNU (обычно меньше: d=64, 2 слоя),
    но на выходе одно число — сколько CX ещё понадобится до цели."""
    arch = "valueu"

    def __init__(self, d=64, de=32, layers=2):
        super().__init__(d=d, de=de, layers=layers)
        self.vhead = mlp(2 * d, d, 1)

    def forward(self, x):
        N = (x.shape[1] // 2).bit_length() - 1
        node, edge = self.features(x, N)
        B, d = x.shape[0], self.d
        h = self.node_in(node)
        e = self.edge_in(edge)
        ce = self.edge_proj(e).view(B, N, N, self.layers_n, d)
        mask = (1 - torch.eye(N, device=x.device)).view(1, N, N, 1)
        for L in range(self.layers_n):
            m = torch.relu(self.a[L](h).unsqueeze(2) + self.b[L](h).unsqueeze(1) + ce[:, :, :, L]) * mask
            h = h + self.upd[L](torch.cat([h, self.out[L](m.sum(2) / max(1, N - 1))], -1))
        return self.vhead(torch.cat([h.mean(1), h.max(1).values], -1)).squeeze(-1)
