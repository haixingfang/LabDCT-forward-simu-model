%dip_binaryrandomvariable   Binary random variable generator.
%    out = dip_binaryrandomvariable(input, p10, p01)
%
%   input
%      Boolean number (1 or 0).
%   p10
%      Real number.
%   p01
%      Real number.

% (C) Copyright 1999-2000               Pattern Recognition Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Cris Luengo, February-May 1999.

%   FUNCTION:
%The binary random variable is generated by altering the input value, if the
%value of a generated random variable is higher than the p10 probability, if
%input is DIP_TRUE, or higher than p01 otherwise.
%ARGUMENTS
%
%  DIPlib      SCIL-Image     Description
%  dip_Random *random            Pointer to a random value structure
%  dip_Boolean input    int input      Input
%  dip_float p10     double p10     Probability of a one to zero transition
%  dip_float p01     double p01     Probability of a one to zero transition
%        int display    Display the return value
%
%EXAMPLE
%Get a binary random variable as follows:
%
%   dip_Random random;
%   dip_float p10, p01, value;
%   p10 = 0.1;
%   p01 = 0.2;
%   DIPXX( dip_BinaryRandomVariable( &random, 1, p10, p01,  &value));
%
%SEE ALSO
% RandomVariable , RandomSeed , UniformRandomVariable , GaussianRandomVariable , PoissonRandomVariable , BinaryRandomVari%able
