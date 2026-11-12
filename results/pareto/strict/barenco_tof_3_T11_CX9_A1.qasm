OPENQASM 3;
include "stdgates.inc";
qubit[6] q;
bit[1] c;
t q[2];
h q[5];
h q[4];
t q[5];
cx q[1], q[5];
tdg q[5];
cx q[0], q[5];
t q[5];
cx q[1], q[5];
tdg q[5];
h q[5];
s q[2];
t q[5];
cx q[2], q[5];
t q[5];
s q[4];
t q[4];
cx q[5], q[4];
s q[4];
t q[4];
cx q[2], q[4];
t q[4];
cx q[5], q[4];
t q[4];
cx q[2], q[4];
h q[4];
h q[5];
c[0] = measure q[5];
if (c[0]) {z q[2];}
if (c[0]) {z q[3];}
if (c[0]) {cz q[0], q[1]; z q[3]; x q[5];}
