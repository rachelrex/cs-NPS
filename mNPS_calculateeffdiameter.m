function [effD] = mNPS_calculateeffdiameter( N, L, d)
% N= average deltaI/I
% L= overall length of the channel
% d= known diameter of polystyrene beads 
% De= variable for effective diameter
% effD= actual calculated effective diameter
syms De
A = (L*N)/(d^3);
B = -1;
C = -0.8*L*N;
eqn= A*(De^3) + B*De + C == 0;

effD=vpasolve(eqn,De);