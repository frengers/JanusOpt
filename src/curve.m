%Slope that is grass lined except for a bare soil zone of 40m
clear all
close all

%Add just the regular file folder
%addpath(genpath('/Users/francisrengers/Documents/MATLAB/BijouCreek/1D Knickpoint&Erosion Model/Sept2012Paper/VaryRain_FluDepGrassOpaquon'));
%cd '/Users/francisrengers/Documents/MATLAB/BijouCreek/1D Knickpoint&Erosion Model/Sept2012Paper/VaryRain_FluDepGrassOpaquon'
%Add the export_fig folder
%addpath(genpath('/Users/francisrengers/Documents/MATLAB/BijouCreek/1D Knickpoint&Erosion Model/Sept2012Paper/VaryRain_FluDepGrassOpaquon/export_fig'));
%cd '/Users/francisrengers/Documents/MATLAB/BijouCreek/1D Knickpoint&Erosion Model/Sept2012Paper/VaryRain_FluDepGrassOpaquon/export_fig'

opt_var_width = 1;  % 0 or 1
width_at_1km2 = 4;  % m

load 'StrdxArray0.txt';
load 'StrElevArray0.txt';
load 'xcell.mat'
in =importdata('InputHack40.txt');

set(0, 'DefaultAxesFontName', 'Arial')
tic;%recording the time

tmax =2;%s
Time=in.data(1);%years
events=10;%number of rainfall events per year
minperevent=30/(60*24*365);%yr The minutes per rainfall event in years
frac_fluvial_erosion_time = events*minperevent; %GT
%GT Tmax=Time*minperevent*events;
Tmax = Time; %GT
dT=in.data(2);%year

for loop =1:1
    if loop==1
        %line='k'
        line=[105/255 105/255 105/255];%This is dark gray
        avgobsq=0.025;%m^3/s This is an average discharge from my flume.
        flumew=(9*2.5)/100;%m 9in*2.5cm/100cm to get flume width
        qin=(avgobsq/flumew);%m/s This is the average discharge measured from last summer
        Agglat=1e-3;%m/yr
        dxsize=1;
        hcheight=4;
        %[newx,newz,zcell,xcell]=IniprofileSameUp_gt(StrdxArray0, StrElevArray0, dxsize, hcheight);
        %xcellsort=sort(xcell,'descend');%m
        zcell = (0.0001192.*(xcell).^2)+(-0.1157.*xcell)+1026;
        zcell(end-70:end) = zcell(end-70:end)-2;
        maxrows=length(2:tmax)*length(1:dT:Tmax)+1;
        maxcols=length(xcell);
        Originalz=zcell;
        shape='-';
        Retreattype='Constant';
        Ptype='Constant';
        Pmmphr=in.data(8);%mm/hr Precip
    elseif loop==2
        Pmmphr=50;%mm/hr Precip
        line=[190/255 190/255 190/255]%This is med lighter gray
        avgobsq=0.025;%m^3/s This is an average discharge from my flume.
        flumew=(9*2.5)/100;%m 9in*2.5cm/100cm to get flume width
        qin=(avgobsq/flumew);%m/s This is the average discharge measured from last summer
        Agglat=1e-3;%m/yr
        dxsize=1;
        hcheight=4;
        [newx,newz,zcell,xcell]=IniprofileSameUp(StrdxArray0, StrElevArray0, dxsize, hcheight);
        xcellsort=sort(xcell,'descend');%m
        maxrows=length(2:tmax)*length(1:dT:Tmax)+1;
        maxcols=length(xcell);
        Originalz=zcell;
        shape='-';
        Retreattype='Constant'
        Ptype='Constant'
    elseif loop==3
        Pmmphr=70;%mm/hr Precip
        line='r'
        %line=[190/255 190/255 190/255]%This is med lighter gray
        avgobsq=0.025;%m^3/s This is an average discharge from my flume.
        flumew=(9*2.5)/100;%m 9in*2.5cm/100cm to get flume width
        qin=(avgobsq/flumew);%m/s This is the average discharge measured from last summer
        Agglat=1e-3;%m/yr
        dxsize=1;
        hcheight=4;
        [newx,newz,zcell,xcell]=IniprofileSameUp(StrdxArray0, StrElevArray0, dxsize, hcheight);
        xcellsort=sort(xcell,'descend');%m
        maxrows=length(2:tmax)*length(1:dT:Tmax)+1;
        maxcols=length(xcell);
        Originalz=zcell;
        shape='-';
        Retreattype='Constant'
    end
    
    %%
    %Constants
    %%
    %Constants
    
    %Pmmphr=50;%Rain in mm/hr
    Immphr=20;%Infiltration in mm/hr
    vs=3.47e-3;%m/s settling velocity of silt
    xstar=1;%m
    %Sediment Discharge Initial Conditions
    qs(1,2)=0;%Sed. discharge coming in the first cell
    qslatper=0;%Percent of the equilibrium sed. conc.(just another way to tweek the equilibrium sed. concentration) If you don't want to use this leave value as 1.
    %Agglat=0;%9e-6;%m/s This represents sed. coming in laterally from the channel sides
    %Critical Shear Stess
    tauc=in.data(3);%kg/(m*s^2)=Pascals from Prosser and Dietrich 1995
    taucWepp=in.data(4);%kg/(m*s^2)=Pascals from WePP Compendium Opequon Soil
    CRetreat=0.5;%m/yr This is the headcut retreat rate
    knickthresh=2.0;%m
    g=9.81;%m/s^2
    rho=1000;%kg/m^3
    sigma=2650;%kg/m^3
    mincon=60;%sec/min change to 1 to convert to time in sec, must change in conjuction with yearconversion
    yrcon=525948.766;%min/year change to 1 to convert to time in sec, must change in conjuction with minconversion
    m=1;% exponent for erosion law (tau-tauc)^m
    phi=0.5;%Porosity
    k=3.37e-3;%s/m WEPP Compendium
    
    
    tauc = 280; %GT FOR DEBUG!

    
    %%
    %Setting up a profile
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Setting up the Profile
    %    dxsize=1;
    %    [newx,newz,zcell,xcell,Getupslp]=Iniprofile(StrdxArray0, StrElevArray0,  dxsize);
    %    xcell=xcell;
    %xcell=0.1:0.1:10;
    %xcell=1:10;
    %     Slope=3e-2;
    %     zupall(1,1)=25;
    %     for i=2:5
    %         zupall(1,i)=zupall(1,i-1)-Slope*(1);
    %     end
    %     zdownall(1,1)=zupall(1,1)-5;
    %     for i=2:5
    %         zdownall(1,i)=zdownall(1,i-1)-Slope*(1);
    %     end
    %     zcell=[zupall zdownall];
    dx=[(xcell(1,2)-xcell(1,1)) diff(xcell(1,:))];
    %dx=ones(1,length(xcell));%m
    %     xcell=xcell;%+1000;%m Influences my erosion rate
    %     newx=newx;%+1000;
    %     xcellsort=sort(xcell,'descend');%m
    maxrows=length(2:tmax)*length(dT:dT:Tmax)+1;
    maxcols=length(xcell);
    
    
    
    %%
    %Initial Conditions
    
    %%
    %Initial Conditions
    
    %Precipitation
    Punits=((Pmmphr/1000)/3600);%Rain in m/s
    I=((Immphr/1000)/3600);%Infiltration in m/s
    
    %Geometery: Channel Width, Slope, Picking Headedge/Headcell
    B=dx(1,1);%(m) channel width. Unless I'm creative about this, I can't easily use a hydraulic geomometry relationship (e.g. B=Q^0.5) as this becomes circular.
    dz0=[(zcell(1,2)-zcell(1,1)) diff(zcell(1,:))];%By taking the difference of the elevation values over the dx, I am using the upwind differencing method
    headedge =find(dz0<=-knickthresh);
    %headcell is the cell at the lip of the headcut
    headcell=headedge-1;
    Stop=-dz0(1,1:headedge-1)./dx(1,1:headedge-1);%Slope above headcut
    Sbtm=-dz0(1,headedge+1:end)./dx(1,headedge+1:end);%Slope below headcut
    zcell0(1,:)=zcell(1,:);%m
    xedge=xcell-dx./2;%m
    xedgeOri=xedge;
    
    %Initial Plot &Plotting parameter
    
    %Turn this on to look add the long profile of Bijou to the plot
    %plot(newx,newz,'b--')
    %hold on
    
    i=1;
    l=1;
    counter=1;
    count=1;
    
    %constq=ones(1,length(xcell));
    %q=qin.*constq;
    cumero=0;%This is the cumulative measure of R*dT
    %just initiallizing some values
    %dzdt=zeros(maxrows,maxcols-1);
    %%GT qs=zeros(maxrows,maxcols-1);
    qs=zeros(1,maxcols-1);  %%GT
    %tauavg=zeros(maxrows,maxcols-1);
    
    % GT moved these out of the loop, to improve performance (only needs to be done once)
    hackc = 5.23;%This is the coefficient I'm using for Hacks law
    hexp = 0.4286;%This is the hacks coefficient
    hcthresh=0.1; %GT: why not 0.1 m?
    
    % GT moved this too: no need to recalculate every time step unless using
    % stochastic or time-varying ppt
    drainage_area = (1/hackc)*xedge.^(1/hexp);
    % GT added option for variable width
    if opt_var_width
        B = width_at_1km2*sqrt(drainage_area./1e6);
    else
        B = B*ones(size(q));
    end
    q=((Punits-I)./B).*drainage_area;%Hacks: Length=aA^h where a is constant, A is drainage area, and h is 0.5-0.6.
    R = CRetreat;  %GT: assume it's R unless told otherwise
    %GT curious about plotting the minimum height at each cell during the
    %run: this is bedrock surface
    bedrock_surface = zcell;
    
     %GT: let's try height dependent
    Retreattype = 'HeightDept';
    if Retreattype=='HeightDept'
        hc_retreat_rate = zeros( 1, ceil( Tmax/dT ) );
        time_step = 1;
    end
    
    %%
    %Huge for loop
    %%
    [r tslice]=size(dT:dT:Tmax);
    stillretreat=1;
    savetime = dT;
    for T=dT:dT:Tmax
        %if T==10+dT;
        %             keyboard
        %         end
        switch Ptype
            %GT case 'Constant'
                %GTP=Punits;
            case 'Stochastic'
                P=exprnd(Punits);
                % GT moved a copy of the line below into here:
                q=((P-I)/B).*(1/hackc)*xedge.^((1/hexp));%Hacks: Length=aA^h where a is constant, A is drainage area, and h is 0.5-0.6.
                if single(T)==single(savetime);
                    psave(count,1)=P;
                    count=count+1;
                end
        end
        
        switch Retreattype
            %GT case 'Constant'
            %GT     R=CRetreat;
            case 'DischargeDept'
                %This is the discharge dependent retreat rate
                qmax=max(q);
                qfrac=q./qmax;
                R=qfrac(headcell)*2*CRetreat;
            case 'HeightDept'
                %This is the height dependent retreat rate
                hmax=6;%m
                cur_ht=zcell(headcell)-zcell(headcell+1);
                hfrac=cur_ht/hmax;
                R=hfrac*2*CRetreat;
                headcut_retreat_rate(time_step) = R;
                time_step = time_step + 1;
        end
        
        %     HEADCUT LOOP
        if ((zcell(headcell)-zcell(headcell+1)>hcthresh) && headedge>3) %If the headcut is > 10cm continue to retreat, else stop retreating
            %GT [headedge headcell xcell zcellknick xedge cumero]=knickpointfunvDepSlope(xcell,zcell,i,(Time/tslice),R,knickthresh,xedge,headedge,headcell,xstar,Stop,dx,xedgeOri,cumero, Sbtm, hcthresh,stillretreat);
            [headedge headcell xcell zcellknick xedge cumero]=knickpointfunvDepSlope(xcell,zcell,i,dT,R,knickthresh,xedge,headedge,headcell,xstar,Stop,dx,xedgeOri,cumero, Sbtm, hcthresh,stillretreat);
            zcell=zcellknick(1,:);
        %elseif (zcell(headcell)-zcell(headcell+1)<=hcthresh) && (stillretreat==1)
            %[headedge headcell xcell zcellknick xedge cumero]=knickpointfunvDepSlope(xcell,zcell,i,(Time/tslice),R,knickthresh,xedge,headedge,headcell,xstar,Stop,dx,xedgeOri,cumero, Sbtm, hcthresh,stillretreat);
        else
            zcell=zcellknick(1,:);
            stillretreat=0;
        end
        
        dz=[(zcell(1,2)-zcell(1,1)) diff(zcell(1,:))];
        dx=[(xcell(1,2)-xcell(1,1)) diff(xcell(1,:))];
        epsilon=1e-2;
        Stop=-dz(1,1:headedge-1)./dx(1,1:headedge-1);%Slope above headcut
        %Stop(Stop<0)=epsilon;
        Sbtm=-dz(1,headedge+1:end)./dx(1,headedge+1:end);%Slope below headcut
        %Sbtm(Sbtm<0)=epsilon;
        S=[Stop Sbtm];
        i=i+1;%Here i is a proxy for time. It records how many times we go through the loop
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Water
        %Note when q gets too big, you go unstable
        %IF YOU USE H WITH CHANGING SLOPE, MAKE HTOP=0 OR MAKE SURE S(1) IS NOT
        %0
        lenzone=in.data(5);%m downstream of headcut of material that is different
        n=in.data(6);%s/m^3 This is manning's n from Chow 1ta959 Pasture High Grass Maximum
        nbare=in.data(7);%s/m^3 This is manning's n from Chow 1959 Earth clean after weathering
        htop=((q(1,1:headedge-1)*n)./(abs(Stop(1:end)).^(1/2))).^(3/5);%m
        hbtm1=((q(1,headedge+1:headedge+lenzone)*nbare)./(Sbtm(1:lenzone).^(1/2))).^(3/5);%m
        hbtm2=((q(1,headedge+lenzone+1:end)*n)./(Sbtm(lenzone+1:end).^(1/2))).^(3/5);%m
        hbtm=[hbtm1 hbtm2];%m
        htopcell=[((diff(htop)./2)+htop(1:end-1)) htop(end)];%m Water Height in cells (avg. water height from edges)
        hbtmcell=[hbtm(1) ((diff(hbtm)./2)+hbtm(1:end-1))];%m Water Height in cells (avg. water height from edges)
        hcell=[htopcell hbtmcell];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Erosion
        %%%Erosion at the edges
        tautop=rho*g.*htop.*Stop;%kg/(m*s^2)
        taubtm=rho*g.*hbtm.*Sbtm;%kg/(m*s^2)
        
        %%%%Erosion in the cells
        tautopavg=[((diff(tautop)./2)+tautop(1:end-1)) tautop(end)];%kg/(m*s^2) Shear is calculated on the edge. So to get the shear in the cell I just take the average shear of the edges. I am making the shear in the last cell equal to the shear on the last edge.
        taubtmavg=[taubtm(1) ((diff(taubtm)./2)+taubtm(1:end-1))];%kg/(m*s^2) Average edge shear stress per each cell.
        tauavg=[tautopavg taubtmavg];%Recording purposes only
        Ecelltopex=(tautopavg-tauc);
        Ecelltopex(Ecelltopex<0)=0;%Get rid of any excess shear that is NEGATIVE
        Ecelltop=(k.*(Ecelltopex).^m).*mincon.*yrcon.*(1/sigma);%[m/yr] Erosion
        
        Ecellbtmex1=(taubtm(1:lenzone)-taucWepp);%kg/(m*s^2)
        Ecellbtmex2=(taubtm(lenzone+1:end)-tauc);%kg/(m*s^2)
        Ecellbtmex=[Ecellbtmex1 Ecellbtmex2];
        Ecellbtmex(Ecellbtmex<0)=0;%Get rid of any excess shear that is NEGATIVE
        Ecellbtm=(k.*(Ecellbtmex).^m).*mincon.*yrcon.*(1/sigma);%[m/yr]
        Ecell=[Ecelltop Ecellbtm];%[m/yr]
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Actual Erosion & Deposition, Perturbing the equilibrium condition
        
        %%GT qs(i,1)=0;%[m^3/yr]
        qs(1,1)=0;%[m^3/yr] %%GT
        %%GT Cedge(1,1)=(qs(i,1))/(q(1,1)*mincon*yrcon)*B;%[unitless]
        Cedge(1,1)=(qs(1,1))/(q(1,1)*mincon*yrcon)*B(1);%[unitless] %%GT
        Dcell(1,1)=vs*Cedge(1,1)*mincon*yrcon;%[m/yr]
        dzdt=(Dcell(1,1)-Ecell(1,1))/(1-phi)+Agglat;
        dzdt = dzdt * frac_fluvial_erosion_time;  %GT
        %Now just run for loop for the cells above headcut
        for j=2:(length(xcell)-1)
            %%GT Cedge(1,j)=(qs(i,j-1))/(q(1,j)*mincon*yrcon*B);%[unitless]
            Cedge(1,j)=(qs(1,j-1))/(q(1,j)*mincon*yrcon*B(j));%[unitless]  %%GT
            Dcell(1,j)=vs*Cedge(1,j)*mincon*yrcon;%[m/yr]
            %%GT qs(i,j)=qs(i,j-1)-(Dcell(1,j)*dx(1,j)*B)+(Ecell(1,j)*dx(1,j)*B);%[m^3/yr]
            qs(1,j)=qs(1,j-1)-(Dcell(1,j)*dx(1,j)*B(j))+(Ecell(1,j)*dx(1,j)*B(j));%[m^3/yr] %%GT
            %Here's the deal. if you have tons of erosion, then you
            %will get tons of deposition in the next cell, then your
            %deposition will be huge compared to your erosion and you
            %will get a negative sediment discharge (qs) which is
            %impossible so I reset it to 0 below. Beware that then
            %every other cell will be 0 as you alternate between
            %erosion and deposition.
            %%GT if qs(i,j)<1e-5
            %%GT    qs(i,j)=0;
            if qs(1,j)<1e-5  %%GT  BE CAREFUL OF STUFF LIKE THIS, IT OFTEN BITES IN THE NETHERS!
                qs(1,j)=0;   %%GT
            end
            dzdt(1,j)=(Dcell(1,j)-Ecell(1,j))/(1-phi)+Agglat;
            dzdt(1,j)=dzdt(1,j)*frac_fluvial_erosion_time;  %GT
            %depratio(i,j)=Dcell(i,j)/Ecell(i,j);%For Debugging. Ratio of Dep. to Erosion
        end
        
        zcell(1,1:end-1)=zcell(1,1:end-1)+dzdt*dT;
        dz=[(zcell(1,2)-zcell(1,1)) diff(zcell(1,:))];
        Stop=-dz(1,1:headedge-1)./dx(1,1:headedge-1);%Slope above headcut
        %Stop(Stop<0)=epsilon;
        Sbtm=-dz(1,headedge+1:end)./dx(1,headedge+1:end);%Slope below headcut
        %Sbtm(Sbtm<0)=epsilon;
        S=[Stop Sbtm];
        zchange=zcell-Originalz;
        %if rem(Tmax,10*dT)==0
        
        % Keep track of bedrock surface as lowest of erosion found so far
        bedrock_surface( zcell < bedrock_surface ) = zcell( zcell < bedrock_surface );
        
%         if T>=next_plot  % GT
%             plot(xcell,zcell)
%             hold on
%             softspotleft = headcell;
%             softspotright = min( headcell+lenzone, length(xcell) );
%             plot(xcell(softspotleft:softspotright),zcell(softspotleft:softspotright),'r');
%             plot(xedge(1:end-1),(qs/500+1000),'c')
%             plot(xcell,bedrock_surface,'k')
%             hold off
%             drawnow
%             next_plot = next_plot + plot_interval;
%         end
        
  
        if single(T)==single(savetime);
           savetime = single(T+(in.data(9)*dT));
           %This is where I need to save the data
           zcellsave(counter,:)=zcell;
           savek(counter,:)=k;
           counter=counter+1;
        end
        
    end
end

save('zcell.mat','zcellsave')
lastcell = zcell;
save('lastz.mat','lastcell')
save('xcell.mat','xcell')

%%
elapsedtime=toc

%Required for Torque Only
exit
