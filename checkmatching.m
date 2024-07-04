x=7;
matching=pulses(1:x,1:8)==prepulses(1:x,1:8)
close=abs(pulses(1:x,1:8)-prepulses(1:x,1:8))<=2