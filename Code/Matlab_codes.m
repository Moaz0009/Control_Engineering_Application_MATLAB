 % matlap engineering control codes
 % ===================================================================
 %laplace transform of impulse input
 syms  t
 x=laplace(dirac(t))
 %===================================================
 %laplace transform of unit step input
 syms  t
 x=laplace(5*t^0)
 %===================================================
 %laplace transform of unit ramp input
 syms  t
 x=laplace(5*t)%ramp function
 %===================================================
 %laplace transform of sin input
 syms  t
 x=laplace(sin(5*t)) %sin function
 %===================================================
 %laplace transform of exp input
 syms  t
 x=laplace(exp(-5*t))%exp function
 %===================================================
 % ilaplace Inverse Laplace transform.
 f=ilaplace(x)
 %===================================================
 % closed loop analysis
 % ==================================================
 %  partial fraction from tf
 % (r(0)/s-p(0))+(r(1)/s-p(1))+......+k
 num=[2];
 den=[1 2 2];
 [r,p,k]=residue(num,den)
 % tf from partial fraction
 [num,den]=residue(r,p,k)
 % ===================================================
 % block diagram 
 G1=tf([1 0],[1 1]);
 G2=tf([1],[1 1 1]);
 Hs=tf([1],[1 0]);
 Gs=series(G1,G2);
 % feedback(Gs,Hs)=>negative feed back
 % feedback(Gs,Hs,+1)=>positive feed back
 T_F_c_l=feedback(Gs,Hs,+1)
 %===================================================
 % response (closed loop tf)
 % =================================================
 num=[2];
 den=[1 2 2];
 figure(1)
 % step responce
 step(num,den)
 figure(2)
 %  step responce in 10 sec
 step(num,den,10)
 figure(3)
 %  impulse responce in 10 sec
 impulse(num,den,10)
 %===============================================
 % response to arbitrary input
 num=[2];
 den=[1 2 2];
 t=[0:.1:10];
 % response to ramp input
 inp=t;   %ramp
 figure(1)
 lsim(num,den,inp,t) 
 % response to sinsudal input
 inp=sin(t);   %sin
 figure(2)
 lsim(num,den,inp,t)  
 %================================================
 %tf2zp  Transfer function to zero-pole conversion.
 num=[2];
 den=[1 2 2];
 [z,p,k]=tf2zp(num,den)
 %===============================================
 %zp2tf zero-pole to tf conversion
 [num,den]=zp2tf(z,p,k)
 %===============================================
 % finding poles&zeros
 num=[2];
 den=[1 2 2];
 sys=tf(num,den)
 P= pole(sys) %find  poles 
 z=zero(sys)  %find zeros
 pzmap(sys)  % draw of poles & zeros0
 %==============================================
 % tf to state space
 num=[1];%without gain
 den=[1 2 2];
 [a,b,c,d]=tf2ss(num,den)
 % state space to tf
 [num,den]=ss2tf(a,b,c,d)
 % ===================================================
 % routh's test
 % char.eq coefficients
 coef=[1 2 2];
 % if all roots(Only The Real Part) are The Same Sign Then 
 % The System is Staple
 roots(coef)
 % ===================================================
 % open loop analysis
 % ===================================================
% open loop tf
G1=tf([1 0],[1 1]);
G2=tf([1],[1 1 1]);
Hs=tf([1],[1 0]);
Gs=series(G1,G2);
T_F_O_L=series(Gs,Hs)
% ========================================
num=[1];
den=[1 2 2 1];
% root locus using num&den
figure
rlocus(num,den)
% root locus using state space eq
% state space for open loop tf
[a,b,c,d]=tf2ss(num,den)
figure
rlocus(a,b,c,d)
% putting zeta and Wn in the graph
% sgrid([zeta1,zeta2, ...],[wn1,wn2, ....])
wn=sqrt(2);
zeta1=1/sqrt(2);% closed loop zeta
zeta2=0.337;% used in root locus calculations
sgrid([zeta1,zeta2],wn)
%====================================================
% polar(nyquist) plot
% using num&den
num=[1];
den=[1 2 2 1];
figure
nyquist(num,den)
% using re&img&w
w=0.1:0.01:100;
[re,img,w]=nyquist(num,den,w);
figure
plot(re,img)
% adding a circle with radius one to calculate Pm&Gm%check stability 
hold on
rectangle('Position',pos,'Curvature',[1 1])
%// radius
r = 1;
%// center
c = [0 0];
pos = [c-r 2*r 2*r];
rectangle('Position',pos,'Curvature',[1 1])
axis equal
hold off
% =========================================
% bode plot
% without Gm&Pm
figure
bode(num,den)
% with Gm&Pm
figure
margin(num,den)
% calculating Gm,Pm,Wp,Wg(without plot)
[gm,pm,wp,wg]=margin(num,den)
% note that gm is not calculated in db here
% so real Gm=20log10(gm)
% ======================================
% nichols plot
figure
nichols(num,den)
