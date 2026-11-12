OPENQASM 3.0;
include "stdgates.inc";
qubit[47] q;
bit[42] c;
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
h q[25];
h q[26];
h q[27];
h q[28];
h q[29];
h q[30];
h q[31];
h q[32];
h q[33];
h q[34];
h q[35];
h q[36];
h q[37];
h q[38];
h q[39];
h q[40];
h q[41];
h q[42];
h q[43];
h q[44];
h q[45];
h q[46];
t q[1];
t q[7];
t q[12];
t q[14];
t q[15];
t q[16];
t q[17];
t q[18];
t q[19];
t q[21];
t q[22];
t q[24];
t q[26];
t q[27];
t q[28];
t q[30];
t q[31];
t q[33];
t q[34];
t q[35];
t q[36];
t q[37];
t q[38];
t q[23];
t q[43];
t q[20];
t q[46];
cx q[44], q[42];
t q[42];
cx q[8], q[6];
t q[6];
cx q[20], q[25];
t q[25];
cx q[37], q[13];
t q[13];
cx q[42], q[3];
t q[3];
cx q[27], q[29];
t q[29];
cx q[5], q[1];
t q[1];
cx q[1], q[32];
t q[32];
cx q[5], q[6];
t q[6];
cx q[2], q[6];
t q[6];
cx q[5], q[6];
t q[6];
cx q[41], q[2];
t q[2];
cx q[5], q[9];
cx q[11], q[9];
t q[9];
cx q[41], q[3];
cx q[42], q[3];
t q[3];
cx q[41], q[2];
cx q[2], q[5];
t q[5];
cx q[41], q[46];
cx q[3], q[46];
t q[46];
cx q[42], q[11];
cx q[11], q[9];
cx q[41], q[9];
t q[9];
cx q[10], q[9];
t q[9];
cx q[41], q[9];
cx q[3], q[9];
t q[9];
cx q[10], q[9];
t q[9];
cx q[38], q[11];
cx q[10], q[11];
cx q[41], q[11];
t q[11];
cx q[39], q[11];
t q[11];
cx q[41], q[11];
cx q[3], q[11];
t q[11];
cx q[39], q[11];
t q[11];
cx q[10], q[39];
cx q[42], q[39];
cx q[41], q[39];
t q[39];
cx q[41], q[3];
cx q[3], q[39];
t q[39];
cx q[29], q[30];
cx q[30], q[29];
cx q[29], q[30];
cx q[18], q[19];
cx q[19], q[18];
cx q[18], q[19];
cx q[8], q[9];
cx q[9], q[6];
cx q[9], q[8];
cx q[2], q[9];
cx q[8], q[9];
cx q[15], q[16];
cx q[16], q[15];
cx q[15], q[16];
cx q[8], q[12];
cx q[12], q[6];
cx q[12], q[8];
cx q[8], q[12];
cx q[10], q[12];
cx q[11], q[12];
cx q[38], q[12];
cx q[11], q[10];
cx q[11], q[39];
cx q[10], q[11];
cx q[29], q[31];
cx q[31], q[29];
cx q[29], q[31];
cx q[4], q[7];
cx q[7], q[4];
cx q[7], q[6];
cx q[6], q[7];
cx q[8], q[7];
cx q[29], q[32];
cx q[32], q[29];
cx q[29], q[32];
cx q[28], q[29];
cx q[29], q[28];
cx q[28], q[29];
cx q[14], q[15];
cx q[15], q[14];
cx q[14], q[15];
cx q[23], q[24];
cx q[24], q[23];
cx q[23], q[24];
cx q[37], q[38];
cx q[38], q[37];
cx q[37], q[38];
cx q[35], q[36];
cx q[36], q[35];
cx q[35], q[36];
cx q[1], q[41];
cx q[41], q[1];
cx q[1], q[41];
cx q[2], q[41];
cx q[5], q[41];
cx q[44], q[45];
cx q[45], q[44];
cx q[3], q[45];
cx q[44], q[45];
cx q[43], q[44];
cx q[44], q[43];
cx q[43], q[44];
cx q[39], q[42];
cx q[42], q[10];
cx q[42], q[39];
cx q[1], q[42];
cx q[3], q[42];
cx q[10], q[42];
cx q[39], q[42];
cx q[40], q[42];
cx q[33], q[34];
cx q[34], q[33];
cx q[33], q[35];
cx q[35], q[33];
cx q[33], q[35];
cx q[25], q[26];
cx q[26], q[25];
cx q[25], q[26];
cx q[37], q[39];
cx q[39], q[37];
cx q[37], q[39];
cx q[40], q[39];
cx q[46], q[2];
cx q[2], q[46];
cx q[3], q[2];
cx q[46], q[2];
cx q[3], q[10];
cx q[27], q[28];
cx q[28], q[27];
cx q[27], q[28];
cx q[0], q[5];
cx q[5], q[0];
cx q[0], q[5];
cx q[20], q[21];
cx q[21], q[20];
cx q[4], q[6];
cx q[8], q[6];
cx q[4], q[8];
cx q[8], q[4];
cx q[4], q[8];
cx q[22], q[23];
cx q[23], q[22];
cx q[22], q[23];
cx q[14], q[17];
cx q[17], q[14];
cx q[14], q[17];
cx q[40], q[0];
cx q[0], q[40];
cx q[40], q[0];
cx q[43], q[4];
cx q[4], q[43];
cx q[43], q[4];
cx q[10], q[40];
cx q[40], q[10];
cx q[1], q[40];
cx q[10], q[40];
cx q[43], q[13];
cx q[13], q[43];
cx q[43], q[13];
cx q[18], q[20];
cx q[20], q[18];
cx q[18], q[20];
cx q[18], q[22];
cx q[22], q[18];
cx q[18], q[22];
cx q[43], q[14];
cx q[14], q[43];
cx q[43], q[14];
cx q[33], q[37];
cx q[37], q[33];
cx q[33], q[37];
cx q[46], q[10];
cx q[25], q[27];
cx q[27], q[25];
cx q[25], q[27];
cx q[43], q[18];
cx q[18], q[43];
cx q[43], q[18];
cx q[43], q[25];
cx q[25], q[43];
cx q[43], q[25];
cx q[33], q[43];
cx q[43], q[33];
cx q[33], q[43];
sdg q[0];
s q[1];
s q[1];
sdg q[2];
s q[9];
s q[9];
s q[10];
s q[10];
sdg q[12];
s q[13];
s q[15];
s q[16];
s q[17];
s q[19];
s q[20];
s q[21];
sdg q[22];
s q[23];
s q[24];
s q[26];
s q[31];
s q[33];
s q[34];
s q[36];
sdg q[37];
sdg q[39];
s q[40];
s q[40];
sdg q[41];
s q[43];
s q[43];
s q[46];
s q[46];
cz q[0], q[3];
cz q[0], q[38];
cz q[0], q[39];
cz q[1], q[41];
cz q[1], q[43];
cz q[1], q[46];
cz q[2], q[3];
cz q[2], q[46];
cz q[3], q[12];
cz q[3], q[39];
cz q[3], q[40];
cz q[3], q[42];
cz q[3], q[43];
cz q[3], q[44];
cz q[4], q[45];
cz q[5], q[10];
cz q[6], q[7];
cz q[6], q[46];
cz q[7], q[8];
cz q[7], q[10];
cz q[8], q[9];
cz q[8], q[46];
cz q[9], q[12];
cz q[9], q[39];
cz q[9], q[42];
cz q[9], q[46];
cz q[10], q[11];
cz q[10], q[41];
cz q[11], q[12];
cz q[12], q[13];
cz q[13], q[14];
cz q[14], q[15];
cz q[15], q[16];
cz q[16], q[17];
cz q[17], q[18];
cz q[18], q[19];
cz q[19], q[20];
cz q[20], q[21];
cz q[20], q[22];
cz q[22], q[23];
cz q[23], q[24];
cz q[24], q[25];
cz q[25], q[26];
cz q[26], q[27];
cz q[27], q[28];
cz q[28], q[29];
cz q[29], q[30];
cz q[30], q[31];
cz q[31], q[32];
cz q[32], q[33];
cz q[33], q[34];
cz q[33], q[35];
cz q[35], q[36];
cz q[36], q[37];
cz q[37], q[38];
cz q[38], q[39];
cz q[39], q[40];
cz q[42], q[43];
cz q[42], q[45];
cz q[43], q[44];
h q[5];
c[0] = measure q[5];
if ((c[0])) { s q[11]; }
if ((c[0])) { s q[11]; }
if ((c[0])) { s q[12]; }
if ((c[0])) { s q[13]; }
if ((c[0])) { s q[13]; }
if ((c[0])) { s q[41]; }
if ((c[0])) { s q[46]; }
if ((c[0])) { s q[46]; }
if ((c[0])) { cz q[3], q[11]; }
if ((c[0])) { cz q[7], q[46]; }
if ((c[0])) { cz q[9], q[46]; }
if ((c[0])) { cz q[10], q[41]; }
h q[6];
c[1] = measure q[6];
if ((c[1])) { s q[7]; }
if ((c[1])) { s q[7]; }
if ((c[1])) { s q[8]; }
if ((c[1])) { s q[8]; }
if ((c[1])) { cz q[10], q[46]; }
if ((c[1]) && (c[0])) { z q[46]; }
h q[41];
c[36] = measure q[41];
if ((c[36])) { s q[3]; }
if ((c[36])) { s q[3]; }
if ((c[36])) { s q[46]; }
if ((c[36])) { cz q[1], q[46]; }
if ((c[36])) { cz q[3], q[43]; }
h q[7];
c[2] = measure q[7];
if ((c[2])) { sdg q[8]; }
if ((c[2])) { s q[9]; }
if ((c[2])) { s q[9]; }
if ((c[2])) { s q[46]; }
if ((c[2])) { s q[46]; }
h q[8];
c[3] = measure q[8];
if ((c[3])) { s q[10]; }
if ((c[3])) { s q[10]; }
if ((c[3])) { s q[12]; }
if ((c[3])) { s q[12]; }
if ((c[3])) { s q[39]; }
if ((c[3])) { s q[39]; }
if ((c[3])) { s q[42]; }
if ((c[3])) { s q[42]; }
if ((c[3])) { cz q[10], q[46]; }
if ((c[3]) && (c[0])) { z q[46]; }
h q[9];
c[4] = measure q[9];
if ((c[4])) { s q[12]; }
if ((c[4]) && (c[0])) { z q[12]; }
if ((c[4])) { s q[13]; }
if ((c[4])) { s q[13]; }
if ((c[4])) { cz q[3], q[11]; }
h q[10];
c[5] = measure q[10];
if ((c[5])) { s q[39]; }
if ((c[5])) { s q[39]; }
if ((c[5])) { s q[42]; }
if ((c[5])) { s q[42]; }
if ((c[5])) { cz q[3], q[12]; }
if ((c[5]) && (c[0] ^ c[4])) { z q[3]; }
if ((c[5])) { cz q[3], q[39]; }
if ((c[5])) { cz q[3], q[42]; }
h q[46];
c[41] = measure q[46];
if ((c[41])) { s q[3]; }
if ((c[41])) { cz q[2], q[3]; }
h q[11];
c[6] = measure q[11];
if ((c[6])) { s q[3]; }
if ((c[6])) { s q[3]; }
if ((c[6])) { s q[12]; }
if ((c[6]) && (c[0] ^ c[4])) { z q[12]; }
if ((c[6])) { s q[13]; }
if ((c[6])) { s q[13]; }
if ((c[6])) { s q[40]; }
if ((c[6])) { s q[40]; }
if ((c[6])) { s q[43]; }
if ((c[6])) { s q[43]; }
if ((c[6])) { s q[45]; }
if ((c[6])) { s q[45]; }
if ((c[6])) { cz q[3], q[40]; }
h q[12];
c[7] = measure q[12];
if ((c[7])) { s q[13]; }
if ((c[7])) { s q[14]; }
if ((c[7])) { s q[14]; }
h q[13];
c[8] = measure q[13];
if ((c[8])) { sdg q[14]; }
if ((c[8])) { s q[15]; }
if ((c[8])) { s q[15]; }
h q[14];
c[9] = measure q[14];
if ((c[9])) { s q[15]; }
if ((c[9])) { s q[16]; }
if ((c[9])) { s q[16]; }
h q[15];
c[10] = measure q[15];
if ((c[10])) { s q[16]; }
if ((c[10])) { s q[17]; }
if ((c[10])) { s q[17]; }
h q[16];
c[11] = measure q[16];
if ((c[11])) { s q[17]; }
if ((c[11])) { s q[18]; }
if ((c[11])) { s q[18]; }
h q[17];
c[12] = measure q[17];
if ((c[12])) { sdg q[18]; }
if ((c[12])) { s q[19]; }
if ((c[12])) { s q[19]; }
h q[18];
c[13] = measure q[18];
if ((c[13])) { s q[19]; }
if ((c[13])) { s q[20]; }
if ((c[13])) { s q[20]; }
h q[19];
c[14] = measure q[19];
if ((c[14])) { s q[20]; }
if ((c[14])) { s q[21]; }
if ((c[14])) { s q[21]; }
if ((c[14])) { s q[22]; }
if ((c[14])) { s q[22]; }
h q[20];
c[15] = measure q[20];
if ((c[15])) { s q[21]; }
if ((c[15])) { sdg q[22]; }
if ((c[15])) { cz q[21], q[22]; }
h q[21];
c[16] = measure q[21];
if ((c[16])) { s q[22]; }
if ((c[16])) { s q[23]; }
if ((c[16])) { s q[23]; }
h q[22];
c[17] = measure q[22];
if ((c[17])) { s q[23]; }
if ((c[17])) { s q[24]; }
if ((c[17])) { s q[24]; }
h q[23];
c[18] = measure q[23];
if ((c[18])) { s q[24]; }
if ((c[18])) { s q[25]; }
if ((c[18])) { s q[25]; }
h q[24];
c[19] = measure q[24];
if ((c[19])) { sdg q[25]; }
if ((c[19])) { s q[26]; }
if ((c[19])) { s q[26]; }
h q[25];
c[20] = measure q[25];
if ((c[20])) { s q[26]; }
if ((c[20])) { s q[27]; }
if ((c[20])) { s q[27]; }
h q[26];
c[21] = measure q[26];
if ((c[21])) { sdg q[27]; }
if ((c[21])) { s q[28]; }
if ((c[21])) { s q[28]; }
h q[27];
c[22] = measure q[27];
if ((c[22])) { sdg q[28]; }
if ((c[22])) { s q[29]; }
if ((c[22])) { s q[29]; }
h q[28];
c[23] = measure q[28];
if ((c[23])) { sdg q[29]; }
if ((c[23])) { s q[30]; }
if ((c[23])) { s q[30]; }
h q[29];
c[24] = measure q[29];
if ((c[24])) { sdg q[30]; }
if ((c[24])) { s q[31]; }
if ((c[24])) { s q[31]; }
h q[30];
c[25] = measure q[30];
if ((c[25])) { s q[31]; }
if ((c[25])) { s q[32]; }
if ((c[25])) { s q[32]; }
h q[31];
c[26] = measure q[31];
if ((c[26])) { sdg q[32]; }
if ((c[26])) { s q[33]; }
if ((c[26])) { s q[33]; }
h q[32];
c[27] = measure q[32];
if ((c[27])) { s q[33]; }
if ((c[27])) { s q[34]; }
if ((c[27])) { s q[34]; }
if ((c[27])) { s q[35]; }
if ((c[27])) { s q[35]; }
h q[33];
c[28] = measure q[33];
if ((c[28])) { s q[34]; }
if ((c[28])) { sdg q[35]; }
if ((c[28])) { cz q[34], q[35]; }
h q[34];
c[29] = measure q[34];
if ((c[29])) { sdg q[35]; }
if ((c[29])) { s q[36]; }
if ((c[29])) { s q[36]; }
h q[35];
c[30] = measure q[35];
if ((c[30])) { s q[36]; }
if ((c[30])) { s q[37]; }
if ((c[30])) { s q[37]; }
h q[36];
c[31] = measure q[36];
if ((c[31])) { s q[37]; }
if ((c[31])) { s q[38]; }
if ((c[31])) { s q[38]; }
h q[37];
c[32] = measure q[37];
if ((c[32])) { s q[0]; }
if ((c[32])) { s q[0]; }
if ((c[32])) { sdg q[38]; }
if ((c[32])) { s q[39]; }
if ((c[32])) { s q[39]; }
h q[38];
c[33] = measure q[38];
if ((c[33])) { s q[0]; }
if ((c[33])) { s q[3]; }
if ((c[33])) { s q[3]; }
if ((c[33])) { s q[39]; }
if ((c[33])) { s q[43]; }
if ((c[33])) { s q[43]; }
if ((c[33])) { s q[45]; }
if ((c[33])) { s q[45]; }
if ((c[33])) { cz q[0], q[39]; }
if ((c[33])) { cz q[3], q[40]; }
h q[39];
c[34] = measure q[39];
if ((c[34])) { s q[0]; }
if ((c[34])) { s q[0]; }
if ((c[34])) { s q[42]; }
if ((c[34])) { s q[42]; }
if ((c[34])) { cz q[0], q[3]; }
if ((c[34])) { cz q[3], q[42]; }
if ((c[34]) && (c[6] ^ c[33])) { z q[3]; }
h q[40];
c[35] = measure q[40];
if ((c[35])) { s q[3]; }
if ((c[35])) { s q[3]; }
if ((c[35])) { s q[43]; }
if ((c[35])) { s q[43]; }
if ((c[35])) { s q[45]; }
if ((c[35])) { s q[45]; }
h q[42];
c[37] = measure q[42];
if ((c[37])) { s q[44]; }
if ((c[37])) { s q[44]; }
if ((c[37])) { cz q[1], q[3]; }
if ((c[37]) && (c[36])) { z q[3]; }
h q[45];
c[40] = measure q[45];
h q[43];
c[38] = measure q[43];
if ((c[38])) { s q[3]; }
if ((c[38])) { s q[3]; }
if ((c[38])) { sdg q[44]; }
h q[44];
c[39] = measure q[44];
if ((c[39])) { s q[4]; }
if ((c[39])) { s q[4]; }
if ((c[39])) { cz q[1], q[3]; }
if ((c[39]) && (c[36])) { z q[3]; }
if ((c[35])) { x q[0]; }
if ((c[36])) { x q[1]; }
if ((c[41])) { x q[2]; }
if ((c[40])) { x q[4]; }
