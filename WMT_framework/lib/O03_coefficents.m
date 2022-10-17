function [A1,A2,B1,B2] = O03_coefficents(Chla)

A1 	= 0.571 + 0.025 .*  log(0.149 .* Chla);
A2	= 0.223 + 0.010 .*  log(2.329 .* Chla);
B1	= 0.015 + 0.176 .* sqrt(0.462 .* Chla);
B2	= 0.688 + 0.060 .*  log(0.125 .* Chla);