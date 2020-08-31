function [dL_t,dT_t,dM_t,c1,c1_dot,gamma,gammasin] = Lift_wing_mov(cs,ys,Alpha0In,EtaS,Cmac,MaxAlphaStallIn,time,dy,ThetaWIn,Beta0,U,f,alpha_b)
global SIUnit n m Density Viscosity W b GammaIn ThetaAIn
Density=1.2250;
GammaIn=35;
ThetaAIn=4.4;
Viscosity=power(10,-5);
y=ys
%
%Viscosity=Kinematic Viscosity,m2/s,ft2/s
%Airfoil parameters (constants,defined by Liebeck LPT 110AAirfoil)
Alpha0=(pi/180)*Alpha0In; %Angle of section's zero lift line,deg
%EtaS Leading edge suction efficiency
%Cmac Moment coeff. about areodynamic centre,to find dMac
MaxAlphaStall=(pi/180)*MaxAlphaStallIn;
%Maximum limit of the areofoil's stall angle,determine flow attached/separated
%Wing parameters
%W=Weight of wing/whole vehicle,N,lb
%U=Freestream velocity/Flight speed,m/s,ft/s
%f=Flapping Frequency, Hz
%b=Wingspan,m,ft
%cs=cord length for various section,m,ft(Note:user input is in inch if USC)
%end
ThetaW=(pi/180)*ThetaWIn; % Mean pitch angle of chord with respective to flapping axis
Gamma=(pi/180)*GammaIn; % Maximum flapping amplitude,deg
ThetaA=(pi/180)*ThetaAIn % Mean Pitch angle of flapping axis w.r.t U,7.5 is in deg
%Beta0 Dynamic Twist,deg/ft
%Calculate t,y,c,AR and k for each section based on geometry of wing and number of sections
 %Width of each section,m,ft
%i=1:2:2*n;
%ys=i*0.5*dy; % Coordinate along semispan for each section,m,ft
t=time;
%y=ys(12);
%for r=0:1:length(time)
%c(r+1)=cs(r+1); % Matrix with m same rows of cs
%end
c=cs;
w=2*pi*f; % Flapping frequency in rad/s
k=w*c/(2*U); % Reduced frequency,rad
SurfaceArea=3.22917;% Assuming each section is rectangular,m2,ft2
AR=4.40; %Aspect Ratio=b2/surface area of 2 wings
%2.COMPUTE PLUNGING & TWISTING MOTION.
%Assumption:Root Flapping kinematics with no spanwise bending,Phase diff.of -90 between %plunge&pitch.
%time steps of m, and n number of sections per wing
ThetaAVal=ThetaA;
GammaVal=Gamma;
Beta0Val=Beta0;
h=-(GammaVal*cos(w*t))*y; % Plunging displacment of LE in flapping direction,m,ft
hDot=w*GammaVal*y.*sin(w*t) ; % 1st Derivative of h, m/s,ft/s
hDotDot=(w^2)*GammaVal*y.*cos(w*t); % 2nd derivative of h, m/s,ft/s2
ThetaBar=ThetaAVal+ThetaW+alpha_b;%+alpha_b % Section's mean pitch angle,constant,rad
dTheta=-Beta0Val*y.*(sin(w*t))*(pi/180);% Dynamically varying pitch angle,Theta-ThetaBar,rad
Theta=dTheta+ThetaBar; % Pitch angle of chord w.r.t U,variable, rad
ThetaDot=-w*Beta0Val*y.*cos(w*t)*(pi/180);%1st DerivativeofTheta,ThetaDot=dThetaDot,rad/s
ThetaDotDot=(w^2)*Beta0Val*y.*sin(w*t)*(pi/180); % 2nd Derivative of Theta, rad/s2
%3.COMPUTE ANGLE OFATTACK AND VELOCITIES
q=Theta-ThetaAVal
PlungeVel=hDot.*cos(q) % Plunging Velocity,m/s, ft/s
PitchVel=0.75*c.*ThetaDot % Pitching velocity with radius of rotation of 3/4 to LE,m/s,ft/s
ForwardVel=U*dTheta % Forward Velocity,consider wing,s motion so AOA=dTheta,m/s,ft/s
PlungeVelDot=hDotDot.*cos(q)-hDot.*sin(q).*ThetaDot;
% 1st derivative of plunging velocity,m/s2,ft/s2
PitchVelDot=0.75*c.*ThetaDotDot; % 1st derivative of Pitching velocity,m/s2,ft/s2
ForwardVelDot=U*ThetaDot; % 1st derivative of Forward Velocity,m/s2,ft/s2
Alpha=(PlungeVel+PitchVel+ForwardVel)/U; %relative aoa at 3/4 chord due to wing's motion,rad
AlphaDot=(PlungeVelDot+PitchVelDot+ForwardVelDot)/U; % 1st derivative of alpha,rad/s
%Compute Alpha' flows relative AOA at 3/4 chord
C1=0.5*AR/(2.32+AR); % Constant used to compute C'(k) by Scherer(1968)
C2=0.181+(0.772/AR); % Constant used to compute C'(k)
z=(k.^2)+(C2^2);
Fk=1-C1*(k.^2)./z; % C'(k)=F'(k)+iG'(k)
Gk=-C1*C2*k./z; % C(k)jones=AR/(2+AR)C'(k)
DownWash=2*(Alpha0+ThetaBar)/(2+AR);
% Downwash=downwashVel/U(Anderson 1991),expression(Kuethe&Chow)
% Downwash due to mean lift produced byAlpha0 and ThetaBar
AlphaPrime=(AR/(2+AR))*(Fk.*Alpha+c.*Gk.*AlphaDot./(2*U.*k))-DownWash; %flows relative AOA at 3/4chord
%V relative flow velocity at 1/4 chord,m/s,ft/s
V=((U*cos(Theta)-hDot.*sin(q)).^2 + (U*(AlphaPrime+ThetaBar)-0.5*c.*ThetaDot).^2).^0.5;
Vat=U*cos(Theta)-hDot.*sin(q); % Tangential component of flow vel on airfoil
Van=hDot.*cos(q)+0.5*c.*ThetaDot+U*sin(Theta); % Normal component of flow vel on airfoil
Va=(Vat.^2+Van.^2).^0.5; % Flow velocity on airfoil due to wing's motion
vDot=U*AlphaDot-0.25*c.*ThetaDotDot;
% Rate of change of midchord normal vel component due to wings motion
% Linearised time derivative of Van
%4.DECIDEWHETHER FLOWIS ATTACHED OR SEPARATED OVER A SECTION
%Assumption:Dynamic stall effects are not considered
%No negative alpha' will occur hence no lower limit
%Separated for the section at time t if element in criterion greater than MaxAlphaStall
Criterion=AlphaPrime+ThetaBar-0.75*(c.*ThetaDot/U);
Attached=(Criterion<=MaxAlphaStall);
Separated=(Criterion>MaxAlphaStall);
%5.COMPUTE NORMAL FORCE AND CHORDWISE FORCE FOR ATTACHED FLOWCONDITION
%Compute Normal Force for Attached Flow
Cn=2*pi*(AlphaPrime+Alpha0+ThetaBar); % Normal force coefficient
dNc=Density*U*V.*Cn.*c*dy/2; % Section's normal force due to circulation
dNa=Density*pi*(c.^2).*vDot*dy/4; % Normal force due to apparent mass effect
dNattach=Attached.*(dNc+dNa); % Section's total normal force with attached flow
%Compute Chordwise Force for Attached Flow
dDcamber=-pi*Alpha0*(AlphaPrime+ThetaBar)*Density*U.*V.*c*dy;
%Chordwise force due to chamber
dTs=EtaS*pi*((AlphaPrime+ThetaBar-0.25*c.*ThetaDot/U).^2)*Density*U.*V.*c*dy; %chordwiseforce due to LE suction
Re=U*c/Viscosity; % Reynolds number
%Cdf = 0.012;
Cdf=0.89./((log10(Re)).^2.58); % Drag coefficient for turbulent boundary layer by Hoerner(1965)
dDf=Cdf.*Density.*(Vat.^2).*c*dy/2; % Chordwise friction drag
dFx=Attached.*(dTs-dDcamber-dDf); % Section's total chordwise force with attached flow
%6.COMPUTE NORMAL FORCE AND CHORDWISE FORCE FOR SEPARATED FLOWCONDITION
%Assuming totally separated flow occurs abruptly, all chordwise forces are negligible
%Compute Normal Force for Separated Flow
Ccfd=1.98; % Crossflow drag coefficient,for high AR flate plate byHoerner
dNcsep=Ccfd*Density*Va.*Van.*c*dy/2; % Normal force due to circulation for separated flow
dNasep=0.5*dNa; % Normal force to to apparent mass effect for separatedflow
dNsep=Separated.*(dNcsep+dNasep); % Secton's total normal force with separated flow
%7.COMPUTE LIFT, THRUST, POWER & PROPULSIVE EFFICIENCY
dN=dNattach+dNsep;
dL=dN.*cos(Theta)+dFx.*sin(Theta); % Lift acting one each section at diiferent time instants
dT=dFx.*cos(Theta)-dN.*sin(Theta); % Trust acting one each section at diiferent time instant
gammat=GammaVal*cos(w*t); % Dihedral angle at an instant in the flapping cycle,rad
Lt=(2*sum((cos(gammat).*dL)'))'; % Instantaneous lift of entire wings
Tt=(2*sum((dT)'))'; % Instantaneous thrust of the entire wings
AveL=(1/(length(time)+1))*sum(Lt); %Average lift of the wing for one flapping cycle,N,lb
AveT=(1/(length(time)+1))*sum(Tt); %Average thrust of the wing for one flapping cycle,N,lb
%Power for attached flow
%dMa=apparent camber and inertia moments
dMa=-Density*pi*(c.^3)*dy.*((1/16)*ThetaDot*U+(1/128)*c.*ThetaDotDot);
dMac=Cmac*Density*(U^2)*SurfaceArea*c/2;
% Section's pitching moment about aerodynamic centre
dPina=dFx.*hDot.*sin(q)+dNattach.*(hDot.*cos(q)+0.25*c.*ThetaDot)+Attached.*dNa.*(0.25*c.*ThetaDot)-Attached.*(dMac+dMa).*ThetaDot;
dPinsep=dNsep.*(hDot.*cos(q)+0.5*c.*ThetaDot);
% Power for separated flow, dMa and dMac ignored
dPin=dPina+dPinsep; % Power absorbed by each section at different time instants
Pint=(2*sum((dPin)'))'; % Instantaneous aerodynamic power absorbed by the whole wing
AvePin=(1/(length(time)+1))*sum(Pint); %Average input power throughout the cycle,Nm/s,ftlb/s
AvePout=AveT*U; %Average output power from the wing,Nm/s,ftlb/s
if(AvePin==0)
AvePropulsiveEff=0; %Average propulsive efficiency, make sure AvePin!=0
else
AvePropulsiveEff=AvePout/AvePin;
end
Lift=AveL; %N , lb
Thrust=AveT; %N , lb
PowIn=AvePin; %Nm/s=W, ft.lb/s=1.356W
Eff=AvePropulsiveEff;
FlapAxisAngle=ThetaAVal/(pi/180); %deg
MaxFlapAmp=GammaVal/(pi/180); %deg
DynamicTwist=Beta0Val; %m/deg, ft/deg
DM=(sum(dMa'))';
DMac=(sum(dMac'))';
dL_t=2*cos(gammat)*dL;
dT_t = 2*dT;
dM_t = DM+DMac;
c1 = -GammaVal*w*sin(w*t);
c1_dot = -w*w*GammaVal*cos(w*t);
gamma =GammaVal*cos(w*t);
gammasin = GammaVal*sin(w*t);
