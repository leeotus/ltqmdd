OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[11];
h qubits[4];
h qubits[4];
ccx qubits[1],qubits[3],qubits[4];
h qubits[4];
h qubits[4];
x qubits[4];
x qubits[5];
h qubits[6];
h qubits[6];
ccx qubits[4],qubits[5],qubits[6];
h qubits[6];
h qubits[6];
x qubits[5];
x qubits[4];
x qubits[6];
h qubits[9];
h qubits[9];
ccx qubits[6],qubits[7],qubits[9];
h qubits[9];
h qubits[9];
x qubits[6];
x qubits[4];
x qubits[5];
h qubits[6];
h qubits[8];
h qubits[6];
ccx qubits[4],qubits[5],qubits[6];
h qubits[6];
h qubits[8];
ccx qubits[10],qubits[9],qubits[8];
h qubits[8];
h qubits[8];
x qubits[5];
x qubits[4];
h qubits[4];
h qubits[2];
h qubits[4];
ccx qubits[1],qubits[3],qubits[4];
h qubits[4];
cx qubits[8],qubits[0];
cx qubits[8],qubits[1];
cx qubits[8],qubits[5];
h qubits[2];
ccx qubits[0],qubits[1],qubits[2];
h qubits[2];
h qubits[2];
x qubits[2];
h qubits[4];
ccx qubits[2],qubits[3],qubits[4];
h qubits[4];
h qubits[4];
x qubits[2];
x qubits[4];
h qubits[6];
ccx qubits[4],qubits[5],qubits[6];
h qubits[6];
h qubits[6];
x qubits[4];
x qubits[6];
h qubits[7];
h qubits[7];
ccx qubits[6],qubits[8],qubits[7];
h qubits[7];
h qubits[7];
x qubits[6];
x qubits[4];
h qubits[6];
h qubits[6];
ccx qubits[4],qubits[5],qubits[6];
h qubits[6];
h qubits[6];
x qubits[4];
h qubits[5];
h qubits[5];
ccx qubits[4],qubits[8],qubits[5];
h qubits[5];
h qubits[5];
x qubits[2];
h qubits[4];
h qubits[4];
ccx qubits[2],qubits[3],qubits[4];
h qubits[4];
h qubits[4];
x qubits[2];
x qubits[2];
h qubits[3];
h qubits[3];
ccx qubits[2],qubits[8],qubits[3];
h qubits[3];
h qubits[3];
x qubits[2];
h qubits[2];
h qubits[2];
ccx qubits[0],qubits[1],qubits[2];
h qubits[2];
h qubits[2];
x qubits[0];
h qubits[1];
h qubits[1];
ccx qubits[0],qubits[8],qubits[1];
h qubits[1];
h qubits[1];
x qubits[0];
h qubits[8];
h qubits[8];
ccx qubits[10],qubits[9],qubits[8];
h qubits[8];
h qubits[8];
