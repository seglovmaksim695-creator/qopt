# Выпуск 0.8.0 — заметки по пакету (24.11.2026)

- Версия 0.8.0 во всех местах номера: `pyproject.toml`, `qoptct/__init__.py`, `qoptct/qopt_core.py`, `CITATION.cff`, `.zenodo.json` (DOI не менялся), README/README_RU, `tests/test_release.py`.
- Новый extra `projectq = ["projectq==0.8.0"]`; extra `qualtran` обслуживает и `qoptct auto`; core-зависимости не менялись.
- package-data добавлено: `foreign/qiskit/data/*.json`, `foreign/projectq/data/*.json`, `foreign/auto/data/*.json`, `dyadic/data/*.json`, `loop/data/*.json` (реестры). В архив исходников (`MANIFEST.in`) добавлены `results/dyadic|auto|projectq|loop`, скрипты `build_*_registry.py`, `dyadic_runs.py`, `projectq_dyadic_runs.py`, `docs/FOREIGN.md` и документы новых частей.
- Реестры перестраиваются скриптами `scripts/build_{qualtran,projectq,auto,dyadic,loop}_registry.py` (у каждого `--check`); тесты сверяют их с файлами результатов (в колесе без исходного дерева такие тесты пропускаются). У Qiskit реестр — файл `data/registry.json` (собран агентом Q1, тест `test_registry_matches_families`).
- Абсолютных путей Мака в пакете нет (проверяет `tests/test_release.py`); исходник Binary_ECC в пакет не входит; кэши — в `$QOPT_CACHE`.
- В `NOTICE`: изменённые файлы Qualtran (Apache-2.0, Copyright Google LLC), Binary_ECC (EPL-2.0, не распространяется), ссылка на работу Amy о правилах dyadic.
- Решение Максима о публикации (24.11.2026): всё выходит одной версией 0.8.0 после сборки и аудита; 0.5–0.7 отдельно не выкладываются; разработчикам Qualtran отдельно не сообщаем. Публикует только координатор: релиз на GitHub → `publish.yml` → PyPI.
- Перед публикацией (не сделано): аудит сборки; тег и релиз v0.8.0; PyPI/Zenodo.

---

# Выпуск 0.7.0 — заметки по пакету (23.11.2026)

- Версия 0.7.0 во всех местах номера: `pyproject.toml`, `qoptct/__init__.py`, `qoptct/qopt_core.py`, `CITATION.cff`, `.zenodo.json` (DOI не менялся), README/README_RU, `tests/test_release.py`.
- Новый extra `qualtran = ["qualtran==0.7.0", "z3-solver>=4.12"]` (был в pyproject от I1); core-зависимости не менялись.
- В пакет добавлены данные (package-data): `foreign/qualtran/patches/{manifest.json,src/*.py,diff/*.diff}`, `foreign/qualtran/ecadd/data/*.json`, `foreign/qualtran/data/*.json` (реестр `registry.json` и файл канареек Калиски n=256). В архив исходников (sdist, `MANIFEST.in`) добавлены `results/qualtran/*.json`, `scripts/build_qualtran_registry.py`, `docs/QUALTRAN*.md`.
- Абсолютных путей Мака в пакете нет (проверяет `tests/test_release.py`); кэши — в `$QOPT_CACHE`.
- В `NOTICE` добавлены изменённые файлы Qualtran (Apache-2.0, Copyright Google LLC).
- Реестр перестраивается командой `python scripts/build_qualtran_registry.py` после любого нового результата; тест `test_registry` проверяет, что он совпадает с файлами в `results/qualtran/`.
- Перед публикацией (не сделано): решение Максима о сообщении разработчикам Qualtran; тег и выпуск; PyPI/Zenodo.

---

# Выпуск 0.3.0 как пакет PyPI `qoptct` — заметки (13.11.2026)

Всё сделано локально в `release_qopt/` (ветка main, локальный коммит). **Ничего не отправлено** на GitHub и
**ничего не опубликовано** на PyPI/TestPyPI.

## Что сделано

1. **Новое имя пакета — `qoptct`.** Имя `qopt` на PyPI занято чужим проектом (qutech/qopt, тоже импортируется
   как `qopt`). Поэтому папка `qopt/` переименована в `qoptct/`, а установка, импорт и команда теперь такие:
   `pip install qoptct`, `import qoptct`, команда `qoptct`. Проект и репозиторий по-прежнему называются «qopt»;
   в README, README_RU, CHANGELOG, `__init__` и справке команды это объяснено. Старой команды `qopt` нет
   (чтобы не было конфликта с чужим пакетом). Переменные окружения `QOPT_*` и папка кэша `~/.cache/qopt`
   не менялись. Все внутренние импорты, загрузки модулей, тесты, скрипты и docs обновлены; исторический
   скрипт `results/hardware/exp.py` оставлен как был (он запускался в старом окружении).
2. **PyTorch стал необязательным** (`pip install "qoptct[net]"`). Нужен только для `--net` и `QoptPass(net=True)`;
   без него команда пишет: *the neural mode needs PyTorch… Install it with: pip install "qoptct[net]"* и
   завершается с кодом 2. Раньше обычный режим всё равно загружал torch (через `block_synth.py`) — теперь torch
   загружается только по требованию; это проверяет тест (даже когда torch установлен).
3. **Зависимости:** numpy ≥ 2.0, qiskit ≥ 2.0, pyzx ≥ 0.10 (обязательные). **Python ≥ 3.11** (код это позволяет;
   проверено на 3.11.18, 3.12.16 и 3.14.6).
4. **Плагин Qiskit** (`qoptct/qiskit_qopt.py`, сообщения и документация на английском):
   * `QoptPass` теперь работает с любой схемой: режет её на непрерывные куски Clifford+T, оптимизирует каждый,
     и заменяет кусок, только если новый **доказанно равен** старому, **строго дешевле** (CX, затем T) и
     использует только разрешённые вентили/пары кубитов (basis gates, coupling map, target). Измерения, барьеры,
     повороты и прочее не трогаются. Глобальная фаза сохраняется точно (куски до 10 кубитов).
   * Плагин стадии оптимизации: `transpile(qc, optimization_method="qoptct")` — сначала обычная оптимизация
     Qiskit выбранного уровня (для базиса Clifford+T — специальная стадия Qiskit), затем `QoptPass`.
   * Дополнительно плагин начальной стадии `init_method="qoptct"`: применяет `QoptPass`, пока Тоффоли ещё не
     разложены. Именно здесь основная экономия (см. таблицу ниже). Оба зарегистрированы через entry points
     `qiskit.transpiler.optimization` и `qiskit.transpiler.init` под именем `qoptct`.
5. **Сборка:** `pyproject.toml` (имя, версия 0.3.0, extra `net`, entry points, package-data), `MANIFEST.in`
   (в архив исходников попадают тесты и разрешённые бенчмарки; hwb6 — никогда).
6. **CI:** `.github/workflows/tests.yml` (pytest на Ubuntu и macOS, Python 3.12, против собранного wheel, без
   torch) и `.github/workflows/publish.yml` (публикация на PyPI через Trusted Publishing при публикации релиза
   на GitHub; environment `pypi`, `id-token: write`; проверка, что тег совпадает с версией). Не запускались.
7. **README / README_RU:** установка через pip, раздел «Qiskit plugin» с примером, бейджи PyPI и тестов
   (бейдж PyPI будет «not found», пока пакет не опубликован — это нормально). CHANGELOG: раздел 0.3.0 с пунктом
   о переименовании. CITATION.cff и .zenodo.json: версия 0.3.0.

## Результаты проверок

| проверка | результат |
|---|---|
| `python -m build` | собраны `qoptct-0.3.0-py3-none-any.whl` (13,7 МБ; 16,3 МБ в распакованном виде) и `qoptct-0.3.0.tar.gz` (13,6 МБ) |
| `twine check --strict dist/*` | PASSED, PASSED |
| содержимое wheel | 4 модели `.pt`, 12 файлов `best_clean`, 22 файла `eng1/ours`, 14 программ портфеля; **нет** hwb6, bench_qc, тестов, `.pyc` |
| pytest (тесты из исходников против установленного wheel, без `-e`) | Python 3.12 без torch: **44 passed**; Python 3.11 без torch: **44 passed**; Python 3.14 с torch: **44 passed** |
| `qoptct optimize bench_qc/tof_3.qc` (режим cx) | T/CX 21/18 → **15/12**, доказано |
| `--mode balanced` | записано 8 T / 10 CX (+2 кубита), независимая проверка F = 1 |
| `--mode t` | 8 T / 10 CX (+2), F = 1 |
| `--mode t --clean` и `--mode balanced --clean` | 8 T / 9 CX (+1), F = 1 (модель «чистые анциллы») |
| `--net` без torch | понятная ошибка, код 2; с torch (3.14): mod5_4 → 7 T / 12 CX, метод `rtof3_esop+net`, доказано |
| `transpile` с плагином, равенство операторов | для tof_3, mod5_4, barenco_tof_3 на уровнях 1–3: `Operator(out) == Operator(qc)` (с фазой) — True |

CX/T после `transpile` в базисе Clifford+T (Qiskit 2.6.2):

| схема | обычный Qiskit | `optimization_method="qoptct"` | `+ init_method="qoptct"` |
|---|---|---|---|
| tof_3 | 18 / 19 | 16 / 15 | 12 / 15 |
| mod5_4 | 28 / 22 | 25 / 16 | 11 / 8 |
| barenco_tof_3 | 24 / 24 | 20 / 16 | 16 / 16 |

С картой связей «линия» плагин никогда не ставит CX между несоседними кубитами (проверено тестом); равенство
проверяется `Operator.from_circuit`.

Ограничения: команда и проход используют сигнал `SIGALRM`, поэтому работают на macOS и Linux, но не на Windows
(в классификаторах указаны только macOS и Linux). Wheel большой (13,7 МБ) из-за моделей — это в пределах
лимита PyPI (100 МБ).

## Что нужно сделать вам (пошагово, без спешки)

Нужно один раз «разрешить» нашему репозиторию на GitHub выкладывать пакет на PyPI. Пароли и ключи никуда
копировать не надо.

1. **Регистрация на PyPI.** Откройте <https://pypi.org/account/register/>, заполните имя, адрес почты, имя
   пользователя и пароль. Подтвердите почту по ссылке из письма.
2. **Двухфакторная защита (обязательна на PyPI).** Войдите на pypi.org → справа вверху ваше имя →
   «Account settings» → раздел «Two factor authentication» → включите вход через приложение-генератор кодов
   (например, Google Authenticator или 1Password) и сохраните резервные коды в надёжном месте.
3. **Разрешение для GitHub (Trusted Publisher).** На pypi.org: ваше имя → «Your projects» → слева
   «Publishing» (или прямая ссылка <https://pypi.org/manage/account/publishing/>). Внизу страницы раздел
   «Add a new pending publisher» → вкладка **GitHub**. Впишите ровно так:
   * PyPI Project Name: `qoptct`
   * Owner: `seglovmaksim695-creator`
   * Repository name: `qopt`
   * Workflow name: `publish.yml`
   * Environment name: `pypi`

   Нажмите «Add». Проект `qoptct` появится на PyPI сам при первой публикации; заранее создавать его не нужно.
   («Pending» значит «ожидающий»: имя закрепится за вами, когда пройдёт первая выкладка.)
4. **Сообщите координатору, что шаги 1–3 сделаны.** Больше ничего делать не нужно.

Что сделает координатор после этого:
* отправит локальный коммит на GitHub (push);
* на GitHub в репозитории появится «окружение» `pypi` (его можно создать в Settings → Environments; если его
  нет, GitHub создаст его автоматически при первом запуске);
* создаст и опубликует релиз **v0.3.0** на GitHub. Это запустит `publish.yml`: GitHub сам соберёт пакет,
  проверит, что версия совпадает с тегом, и выложит его на PyPI. Через пару минут пакет будет доступен:
  `pip install qoptct`, страница <https://pypi.org/project/qoptct/>, и бейдж PyPI в README оживёт.
* Если на GitHub подключена интеграция с Zenodo, релиз также создаст новую версию записи Zenodo (из `.zenodo.json`).

Если что-то пойдёт не так (например, ошибка «invalid-publisher»), почти всегда причина — опечатка в одном из
пяти полей шага 3; её можно исправить на той же странице PyPI и перезапустить задание на GitHub.
