OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[30];
h qubits[4];
h qubits[5];
h qubits[10];
h qubits[11];
h qubits[12];
h qubits[13];
h qubits[18];
h qubits[19];
h qubits[24];
h qubits[25];
h qubits[26];
h qubits[27];
h qubits[28];
h qubits[29];
h qubits[4];
ccx qubits[0],qubits[3],qubits[4];
h qubits[4];
x qubits[1];
h qubits[5];
ccx qubits[1],qubits[2],qubits[5];
h qubits[5];
x qubits[1];
x qubits[0];
h qubits[4];
ccx qubits[0],qubits[2],qubits[4];
h qubits[4];
x qubits[0];
h qubits[5];
ccx qubits[1],qubits[3],qubits[5];
h qubits[5];
h qubits[10];
ccx qubits[6],qubits[9],qubits[10];
h qubits[10];
x qubits[7];
h qubits[11];
ccx qubits[7],qubits[8],qubits[11];
h qubits[11];
x qubits[7];
x qubits[6];
h qubits[10];
ccx qubits[6],qubits[8],qubits[10];
h qubits[10];
x qubits[6];
h qubits[11];
ccx qubits[7],qubits[9],qubits[11];
h qubits[11];
h qubits[18];
ccx qubits[14],qubits[17],qubits[18];
h qubits[18];
x qubits[15];
h qubits[19];
ccx qubits[15],qubits[16],qubits[19];
h qubits[19];
x qubits[15];
x qubits[14];
h qubits[18];
ccx qubits[14],qubits[16],qubits[18];
h qubits[18];
x qubits[14];
h qubits[19];
ccx qubits[15],qubits[17],qubits[19];
h qubits[19];
h qubits[24];
ccx qubits[20],qubits[23],qubits[24];
h qubits[24];
x qubits[21];
h qubits[25];
ccx qubits[21],qubits[22],qubits[25];
h qubits[25];
x qubits[21];
x qubits[20];
h qubits[24];
ccx qubits[20],qubits[22],qubits[24];
h qubits[24];
x qubits[20];
h qubits[25];
ccx qubits[21],qubits[23],qubits[25];
h qubits[25];
h qubits[4];
h qubits[5];
h qubits[10];
h qubits[11];
h qubits[18];
h qubits[19];
h qubits[24];
h qubits[25];
x qubits[4];
h qubits[12];
ccx qubits[4],qubits[10],qubits[12];
h qubits[12];
x qubits[4];
h qubits[13];
ccx qubits[5],qubits[11],qubits[13];
h qubits[13];
x qubits[5];
h qubits[13];
ccx qubits[5],qubits[10],qubits[13];
h qubits[13];
x qubits[5];
h qubits[12];
ccx qubits[4],qubits[11],qubits[12];
h qubits[12];
x qubits[18];
h qubits[26];
ccx qubits[18],qubits[24],qubits[26];
h qubits[26];
x qubits[18];
h qubits[27];
ccx qubits[19],qubits[25],qubits[27];
h qubits[27];
x qubits[19];
h qubits[27];
ccx qubits[19],qubits[24],qubits[27];
h qubits[27];
x qubits[19];
h qubits[26];
ccx qubits[18],qubits[25],qubits[26];
h qubits[26];
h qubits[12];
h qubits[13];
h qubits[26];
h qubits[27];
x qubits[12];
h qubits[28];
ccx qubits[12],qubits[26],qubits[28];
h qubits[28];
x qubits[12];
h qubits[29];
ccx qubits[13],qubits[27],qubits[29];
h qubits[29];
x qubits[13];
h qubits[29];
ccx qubits[13],qubits[26],qubits[29];
h qubits[29];
x qubits[13];
h qubits[28];
ccx qubits[12],qubits[27],qubits[28];
h qubits[28];
h qubits[28];
h qubits[29];
