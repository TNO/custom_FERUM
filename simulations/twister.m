%TWISTER   Uniformly distributed pseudo-random numbers.
%   R = TWISTER(N) returns an N-by-N matrix containing pseudo-random values
%   drawn from a uniform distribution on the unit interval.  TWISTER(M,N) or
%   TWISTER([M,N]) returns an M-by-N matrix.  TWISTER(M,N,P,...) or
%   TWISTER([M,N,P,...]) generates an M-by-N-by-P-by-... array.  TWISTER with
%   no arguments returns a scalar.  TWISTER(SIZE(A)) returns an array the same
%   size as A.
%
%   TWISTER produces pseudo-random numbers using the Mersenne Twister
%   algorithm by Nishimura and Matsumoto, and is an alternative to the
%   built-in function RAND in MATLAB.  It creates double precision values in
%   the closed interval [0, 1-2^(-53)], and can generate 2^19937 - 1 values
%   before repeating itself.
%
%   The sequence of numbers generated is determined by the internal state of
%   the generator.  Since MATLAB resets the state at start-up, the sequence of
%   numbers generated will be the same in each session unless the state is
%   changed.  Setting the generator to different states leads to unique
%   computations, but does not improve any statistical properties.  Setting
%   the generator to the same fixed state allows computations to be repeated.
%
%   TWISTER('state',J), where J is a scalar integer, initializes the state of
%   the generator.  There is no simple connection between the sequence of
%   random numbers generated from TWISTER('state',J) and TWISTER('state',J+1).
%   TWISTER('state',0) resets the generator to its initial state.  J may also
%   be an array of integers with length less than 625.
%
%   S = TWISTER('state') returns a 625-element vector of UINT32 values
%   containing the current state of the uniform generator.
%
%   TWISTER('state',S), where S is the output of TWISTER('state'), sets the
%   state of the generator to S.
%
%    Examples:
%
%       Three ways to initialize TWISTER differently each time:
%          twister('state',sum(100*clock))
%          twister('state',100*clock)
%          twister('state',2^32*rand(n,1)) % where n < 625.
%
%       Generate 100 values, reset the state, and repeat the sequence:
%          s = twister('state');
%          u1 = twister(100,1);
%          twister('state',s);
%          u2 = twister(100,1); % contains exactly the same values as u1
%
%       Generate uniform values from the interval [a, b]:
%          r = a + (b-a).*twister(100,1);
%
%       Generate integers uniform on the set 1:n:
%          r = 1 + floor(n.*twister(100,1));
%
%       Generate standard normal random values using the inversion method:
%          z = -sqrt(2).*erfcinv(2*twister(100,1));
%
%   Mex file derived from a copyrighted C program by Takuji Nishimura and
%   Makoto Matsumoto.
%
%   Reference: M. Matsumoto and T. Nishimura, "Mersenne Twister: A
%   623-Dimensionally Equidistributed Uniform Pseudo-Random Number Generator",
%   ACM Transactions on Modeling and Computer Simulation, Vol. 8, No. 1,
%   January 1998, pp 3--30.
%
%   See also RAND, RANDN.

%   Note:  Initializing TWISTER to the scalar integer state 0 actually
%   corresponds to the C call init_genrand(5489).

%   Written by Peter Perkins, The MathWorks, Inc.
%   Revision: 1.0  Date: 2004/12/23
%   This function is not supported by The MathWorks, Inc.
%
%   Requires MATLAB R13.


% ==================== Copyright for Mersenne Twister ====================

%    A C-program for MT19937, with initialization improved 2002/1/26.
%    Coded by Takuji Nishimura and Makoto Matsumoto.
%
%    Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
%    All rights reserved.
%
%    Redistribution and use in source and binary forms, with or without
%    modification, are permitted provided that the following conditions
%    are met:
%
%      1. Redistributions of source code must retain the above copyright
%         notice, this list of conditions and the following disclaimer.
%
%      2. Redistributions in binary form must reproduce the above copyright
%         notice, this list of conditions and the following disclaimer in the
%         documentation and/or other materials provided with the distribution.
%
%      3. The names of its contributors may not be used to endorse or promote
%         products derived from this software without specific prior written
%         permission.
%
%    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%    A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
%    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
%    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
%    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
%    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
%    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
%    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%
%    Any feedback is very welcome.
%    http://www.math.keio.ac.jp/matumoto/emt.html
%    email: matumoto@math.keio.ac.jp
