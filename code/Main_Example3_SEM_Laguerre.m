
%-- Questo programma risolve il problema agli autovalori della linea di
%   trasmissione, lunga pigreco. Se si mettono entrambe condizioni di
%   Neumann o entrambe di Dirichlet, le radici degli autovalori sono
%   multiple di 1: {1, 2, 3,...}

clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- General parameters
Options.Orthogonality=1; % if 1 functions are orthogonal DA IMPLEMENTARE
Options.BCTHreshold=1e-12; % boundary conditions SVD threshold
Options.OCTHreshold=1e-10; % boundary conditions SVD threshold
Options.CCTHreshold=1e-12; % boundary conditions SVD threshold

%-- Patch 1
PatchInfo(1).Ntot=8;
PatchInfo(1).Nquad=57;
PatchInfo(1).L=500;
PatchInfo(1).BCInfo={'Natural','Continuity'};
PatchInfo(1).kFlux=1; % for flux continuity (kFlux=1 means Derivability)
PatchInfo(1).PolyType='Laguerre';
%-- Patch 2
PatchInfo(2).Ntot=6;
PatchInfo(2).Nquad=56;
PatchInfo(2).L=50;
PatchInfo(2).BCInfo={'Continuity','Continuity'};
PatchInfo(2).kFlux=4; % for flux continuity (kFlux=1 means Derivability)
PatchInfo(2).PolyType='Legendre';
%-- Patch 3
PatchInfo(3).Ntot=7;
PatchInfo(3).Nquad=55;
PatchInfo(3).L=500;
PatchInfo(3).BCInfo={'Continuity','Natural'};
PatchInfo(3).kFlux=1; % for flux continuity (kFlux=1 means Derivability)
PatchInfo(3).PolyType='Laguerre';

indMode=1;
k0=2*pi/85;

% [StrFun,Nfun]=f_SynthesizeSEMBasisFunctions(z1,z2,L1,L2,Nf_1,Nf_2,BCInfo);
%--- DA QUI IN POI E' LA FUNCTION DI SINTESI

xRight=0;
for indPatch=1:length(PatchInfo)
    
    Nquad=PatchInfo(indPatch).Nquad;

    if(strcmp(PatchInfo(indPatch).PolyType,'Legendre')==1)
        
        Lpatch=PatchInfo(indPatch).L;
        vn=[0:1:PatchInfo(indPatch).Ntot];
        PatchInfo(indPatch).Ntot=PatchInfo(indPatch).Ntot+1; %-- Ntot+1 poly.
        PatchInfo(indPatch).dudx=2/Lpatch;
        PatchInfo(indPatch).dxdu=Lpatch/2;

        [nodes,weights]=quadad(1,Nquad); % Gauss-Legendre nodes
        [u,I]=sort(nodes); % sorting nodes ...
        wu=weights(I); % ... and weights ...

        PatchInfo(indPatch).x=xRight+Lpatch/2.*u+Lpatch/2;
        PatchInfo(indPatch).wx=Lpatch/2.*wu;

        [f,dfdu,d2fdu2]=f_EvalLegendrePolynomials(vn,u); % internal points
        [f_L,dfdu_L,d2fdu2_L]=f_EvalLegendrePolynomials(vn,-1); % left B.C.
        [f_R,dfdu_R,d2fdu2_R]=f_EvalLegendrePolynomials(vn,1); % right B.C.

    elseif(strcmp(PatchInfo(indPatch).PolyType,'Chebyshev')==1)
        
        Lpatch=PatchInfo(indPatch).L;
        vn=[0:1:PatchInfo(indPatch).Ntot];
        PatchInfo(indPatch).Ntot=PatchInfo(indPatch).Ntot+1; %-- Ntot+1 poly.
        PatchInfo(indPatch).dudx=2/Lpatch;
        PatchInfo(indPatch).dxdu=Lpatch/2;

        [nodes,weights]=quadad(1,Nquad); % Gauss-Legendre nodes
        [u,I]=sort(nodes); % sorting nodes ...
        wu=weights(I); % ... and weights ...

        PatchInfo(indPatch).x=xRight+Lpatch/2.*u+Lpatch/2;
        PatchInfo(indPatch).wx=Lpatch/2.*wu;
        
        [f,dfdu,d2fdu2]=f_EvalChebyshevPolynomials(vn,u); % internal points
        [f_L,dfdu_L,d2fdu2_L]=f_EvalChebyshevPolynomials(vn,-1); % left B.C.
        [f_R,dfdu_R,d2fdu2_R]=f_EvalChebyshevPolynomials(vn,1); % right B.C.
        
    elseif(strcmp(PatchInfo(indPatch).PolyType,'Laguerre')==1)
        
        Lpatch=0;
        [u,wu]=f_GaussLaguerreQuad(Nquad); % Laguerre nodes
        vn=[0:1:PatchInfo(indPatch).Ntot];
        wu=exp(u).*wu; % by this way the weight is in the basis functions!
        PatchInfo(indPatch).Ntot=PatchInfo(indPatch).Ntot+1; %-- Ntot+1 poly.
        [Ln,dLndu]=f_EvalLaguerrePolynomials(vn,u); 
        [Ln_L,dLndu_L]=f_EvalLaguerrePolynomials(vn,0);
        U=ones(PatchInfo(indPatch).Ntot,1)*u;
        WU=ones(PatchInfo(indPatch).Ntot,1)*wu;
        f=Ln.*exp(-U/2);
        dfdu=(dLndu-0.5.*Ln).*exp(-U/2);
        
        if(indPatch==length(PatchInfo))
                       
            PatchInfo(indPatch).dudx=1;
            PatchInfo(indPatch).dxdu=1;
            
            PatchInfo(indPatch).x=PatchInfo(indPatch).dxdu.*u+xRight;
            PatchInfo(indPatch).wx=PatchInfo(indPatch).dxdu.*wu;
            
            f_L=Ln_L; dfdu_L=(dLndu_L-0.5.*Ln_L); f_R=zeros(size(f_L));
            d2fdu2=zeros(size(dfdu)); dfdu_R=zeros(size(dfdu_L));
            d2fdu2_L=zeros(size(dfdu_L)); d2fdu2_R=zeros(size(dfdu_L));
            
        elseif(indPatch==1)
            
            PatchInfo(indPatch).dudx=-1; % map x -> -u
            PatchInfo(indPatch).dxdu=-1; % map x -> -u
            
            PatchInfo(indPatch).x=PatchInfo(indPatch).dxdu.*u;
            % Integral should remain EQUAL!!!
            PatchInfo(indPatch).wx=-PatchInfo(indPatch).dxdu.*wu;
            
            f_R=Ln_L; dfdu_R=(dLndu_L-0.5.*Ln_L); f_L=zeros(size(f_R));
            d2fdu2=zeros(size(dfdu)); dfdu_L=zeros(size(dfdu_R));
            d2fdu2_L=zeros(size(dfdu_R)); d2fdu2_R=zeros(size(dfdu_R));
            
        else % Laguerre can be used only on first and/or last patch!!!
            error(['Wrong Patch ',num2str(indPatch) ,' polynomial type'])
        end
                
    else error(['Wrong Patch ',num2str(indPatch) ,' polynomial type'])
        
    end
    
    PatchInfo(indPatch).f=f;
    PatchInfo(indPatch).dfdx=dfdu*PatchInfo(indPatch).dudx;
    PatchInfo(indPatch).d2fdx2=d2fdu2*PatchInfo(indPatch).dudx^2;
 
    PatchInfo(indPatch).f_L=f_L;
    PatchInfo(indPatch).dfdx_L=dfdu_L*PatchInfo(indPatch).dudx;
    PatchInfo(indPatch).d2fdx2_L=d2fdu2_L*PatchInfo(indPatch).dudx^2;
    
    PatchInfo(indPatch).f_R=f_R;
    PatchInfo(indPatch).dfdx_R=dfdu_R*PatchInfo(indPatch).dudx;
    PatchInfo(indPatch).d2fdx2_R=d2fdu2_R*PatchInfo(indPatch).dudx^2;
    
    xRight=xRight+Lpatch;
    PatchInfo(indPatch).Nfun=PatchInfo(indPatch).Ntot;
    
    PatchInfo(indPatch).coeff_mn=eye(PatchInfo(indPatch).Nfun);
    
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary and continuity conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(length(PatchInfo)>1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Left boundary condition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BCInfo=PatchInfo(1).BCInfo{1};
    if(not(strcmp(BCInfo,'Natural'))) %-- If Natural, no explicit enforcement
        if(strcmp(BCInfo,'Dirichlet'))
            fBC=[PatchInfo(1).f_L.']; %-- Dirichlet
        elseif(strcmp(BCInfo,'Neumann'))
            fBC=[PatchInfo(1).dfdx_L.']; %-- Neumann
        else
            error('Wrong left boundary condition')
        end
        [U,S,V]=svd(fBC);
        S=diag(S);
        indBC=find(abs(S./S(1))<=Options.BCTHreshold);
        if isempty(indBC)==1
            indBC=length(S)+1;
        else
            indBC=indBC(1);
        end
        PatchInfo(1).coeff_mn=V(:,indBC:end);
        PatchInfo(1).Nfun=size(PatchInfo(1).coeff_mn,2); % upgrade Ntot for B.C.
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Right boundary condition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BCInfo=PatchInfo(end).BCInfo{2};
    if(not(strcmp(BCInfo,'Natural'))) %-- If Natural, no explicit enforcement
        if(strcmp(BCInfo,'Dirichlet'))
            fBC=[PatchInfo(end).f_R.']; %-- Dirichlet
        elseif(strcmp(BCInfo,'Neumann'))
            fBC=[PatchInfo(end).dfdx_R.']; %-- Neumann
        else
            error('Wrong right boundary condition')
        end
        [U,S,V]=svd(fBC);
        S=diag(S);
        indBC=find(abs(S./S(1))<=Options.BCTHreshold);
        if isempty(indBC)==1
            indBC=length(S)+1;
        else
            indBC=indBC(1);
        end
        PatchInfo(end).coeff_mn=V(:,indBC:end);
        PatchInfo(end).Nfun=size(PatchInfo(end).coeff_mn,2); % upgrade Ntot for B.C.
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Orthogonality
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(Options.Orthogonality)
    for indPatch=1:length(PatchInfo)
        
        if(strcmp(PatchInfo(indPatch).PolyType,'Legendre'))
            
            vMass=2./(2.*(0:1:(PatchInfo(indPatch).Ntot-1))+1);
            MPatch=PatchInfo(indPatch).dxdu.*diag(vMass);
            
        elseif(strcmp(PatchInfo(indPatch).PolyType,'Chebyshev'))

            WX=ones(PatchInfo(indPatch).Ntot,1)*PatchInfo(indPatch).wx;
            MPatch=conj(PatchInfo(indPatch).f.*WX)*PatchInfo(indPatch).f.';
            
        elseif(strcmp(PatchInfo(indPatch).PolyType,'Laguerre'))

            WX=ones(PatchInfo(indPatch).Ntot,1)*PatchInfo(indPatch).wx;
            MPatch=conj(PatchInfo(indPatch).f.*WX)*PatchInfo(indPatch).f.';
            
        end
        
        MPatch=PatchInfo(indPatch).coeff_mn'*MPatch*PatchInfo(indPatch).coeff_mn;

        [U,S,V]=svd(MPatch);
        S=diag(S);
        I=find(S>=Options.OCTHreshold);
        Nfun=PatchInfo(indPatch).Nfun;
        U=U(:,I)./sqrt((ones(Nfun,1)*S(I).'));
        PatchInfo(indPatch).coeff_mn=PatchInfo(indPatch).coeff_mn*U;
        PatchInfo(indPatch).Nfun=size(PatchInfo(indPatch).coeff_mn,2);

    end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Continuity conditions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Building unique set of basis functions
    Ntot=0;
    for indPatch=1:length(PatchInfo)
        PatchInfo(indPatch).indCC=Ntot+(1:PatchInfo(indPatch).Nfun);
        Ntot=Ntot+PatchInfo(indPatch).Nfun;
    end
    fCC=[];
    
    for indCC=1:(length(PatchInfo)-1) %-- Loop on the interfaces
        CCInfo_L=PatchInfo(indCC).BCInfo{2};
        CCInfo_R=PatchInfo(indCC+1).BCInfo{1};
        %-- Check if the continuity condition is properly indicated
        if(not(strcmp(CCInfo_L,CCInfo_R)))
            error(['Wrong continuity definitions between patches ',num2str(indCC),' and ',num2str(indCC+1)])
        end
        fCCRow=zeros(1,Ntot);
        if(strcmp(CCInfo_L,'Continuity'))
            fCC_L=PatchInfo(indCC).coeff_mn.'*PatchInfo(indCC).f_R;
            fCC_R=PatchInfo(indCC+1).coeff_mn.'*PatchInfo(indCC+1).f_L;
            
            fCCRow(PatchInfo(indCC).indCC)=fCC_L.';
            fCCRow(PatchInfo(indCC+1).indCC)=-fCC_R.';
            fCC=[fCC;fCCRow];
            
        elseif(strcmp(CCInfo_L,'Derivability'))
            fCCRowFlux=zeros(1,Ntot);
            
            fCC_L=PatchInfo(indCC).coeff_mn.'*PatchInfo(indCC).f_R;
            fCC_R=PatchInfo(indCC+1).coeff_mn.'*PatchInfo(indCC+1).f_L;
            fCCFlux_L=PatchInfo(indCC).coeff_mn.'*PatchInfo(indCC).dfdx_R.*PatchInfo(indCC).kFlux;
            fCCFlux_R=PatchInfo(indCC+1).coeff_mn.'*PatchInfo(indCC+1).dfdx_L.*PatchInfo(indCC+1).kFlux;
            
            fCCRow(PatchInfo(indCC).indCC)=fCC_L.';
            fCCRow(PatchInfo(indCC+1).indCC)=-fCC_R.';
            fCCRowFlux(PatchInfo(indCC).indCC)=fCCFlux_L.';
            fCCRowFlux(PatchInfo(indCC+1).indCC)=-fCCFlux_R.';
            
            fCC=[fCC;fCCRow;fCCRowFlux];
            
        else
            error(['Wrong continuity definitions between patches ',num2str(indCC),' and ',num2str(indCC+1)])
        end
        
    end
    
    [U,S,V]=svd(fCC);
    S=diag(S);
    indCC=find(abs(S./S(1))<=Options.CCTHreshold);
    if isempty(indCC)==1
        indCC=length(S)+1;
    else
        indCC=indCC(1);
    end
    coeff_mn_CC=V(:,indCC:end);
    
    for indPatch=1:length(PatchInfo)
        PatchInfo(indPatch).coeff_mn=PatchInfo(indPatch).coeff_mn*coeff_mn_CC(PatchInfo(indPatch).indCC,:);
        PatchInfo(indPatch).Nfun=size(PatchInfo(indPatch).coeff_mn,2); % upgrade Ntot for C.C.
    end
else %----------------------------------- 1 patch analysis
    BCInfo=PatchInfo(end).BCInfo;
    if(not(strcmp(BCInfo{1},'Natural')) || not(strcmp(BCInfo{2},'Natural')))
    if(strcmp(BCInfo{1},'Dirichlet'))
        fBC=[PatchInfo(end).f_L.']; %-- Dirichlet
    elseif(strcmp(BCInfo{1},'Neumann'))
        fBC=[PatchInfo(end).dfdx_L.']; %-- Neumann
    elseif(not(strcmp(BCInfo{1},'Natural')))
        error('Wrong left boundary condition')
    end
    if(strcmp(BCInfo{2},'Dirichlet'))
        fBC=[fBC;PatchInfo(end).f_R.']; %-- Dirichlet
    elseif(strcmp(BCInfo{2},'Neumann'))
        fBC=[fBC;PatchInfo(end).dfdx_R.']; %-- Neumann
    elseif(not(strcmp(BCInfo{2},'Natural')))
        error('Wrong left boundary condition')
    end
    [U,S,V]=svd(fBC);
    S=diag(S);
    indBC=find(abs(S./S(1))<=Options.BCTHreshold);
    if isempty(indBC)==1
        indBC=length(S)+1;
    else
        indBC=indBC(1);
    end
    PatchInfo(end).coeff_mn=V(:,indBC:end);
    PatchInfo(end).Nfun=size(PatchInfo(end).coeff_mn,2); % upgrade Ntot for B.C.
    end
    
    if(Options.Orthogonality)
        for indPatch=1:length(PatchInfo)
            
            if(strcmp(PatchInfo(indPatch).PolyType,'Legendre'))
                
                vMass=2./(2.*(0:1:(PatchInfo(indPatch).Ntot-1))+1);
                MPatch=PatchInfo(indPatch).dxdu.*diag(vMass);
                
            elseif(strcmp(PatchInfo(indPatch).PolyType,'Chebyshev'))
                
                WX=ones(PatchInfo(indPatch).Ntot,1)*PatchInfo(indPatch).wx;
                MPatch=conj(PatchInfo(indPatch).f.*WX)*PatchInfo(indPatch).f.';
                            
            elseif(strcmp(PatchInfo(indPatch).PolyType,'Laguerre'))
                
                WX=ones(PatchInfo(indPatch).Ntot,1)*PatchInfo(indPatch).wx;
                MPatch=conj(PatchInfo(indPatch).f.*WX)*PatchInfo(indPatch).f.';
                
            end
            
            MPatch=PatchInfo(indPatch).coeff_mn'*MPatch*PatchInfo(indPatch).coeff_mn;
            
            [U,S,V]=svd(MPatch);
            S=diag(S);
            I=find(S>=Options.OCTHreshold);
            Nfun=PatchInfo(indPatch).Nfun;
            U=U(:,I)./sqrt((ones(Nfun,1)*S(I).'));
            PatchInfo(indPatch).coeff_mn=PatchInfo(indPatch).coeff_mn*U;
            PatchInfo(indPatch).Nfun=size(PatchInfo(indPatch).coeff_mn,2);
            
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test problem: second order transmission line equation (eigvl. problem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=zeros(PatchInfo(1).Nfun); 
K=M; Mtrue=M;

[D1,D2]=f_EvalAnalyticLegendreIntegrals(2*PatchInfo(1).Nfun);

for indPatch=1:length(PatchInfo)
    Weight=ones(PatchInfo(indPatch).Ntot,1)*PatchInfo(indPatch).wx;
    if(strcmp(PatchInfo(indPatch).PolyType,'Chebyshev')) % numerical
        Mi=conj(PatchInfo(indPatch).f.*Weight)*PatchInfo(indPatch).f.';
        Ki=conj(PatchInfo(indPatch).dfdx.*Weight)*PatchInfo(indPatch).dfdx.';
    elseif(strcmp(PatchInfo(indPatch).PolyType,'Legendre')) % analytical
        vn=[0:1:PatchInfo(indPatch).Ntot-1];
        Mi=diag(2./(2.*vn+1))./PatchInfo(indPatch).dudx; % Gram matrix
        Ki=(Mi)*D2(1:PatchInfo(indPatch).Ntot,1:PatchInfo(indPatch).Ntot)*PatchInfo(indPatch).dudx^2;
        BoundaryTerm=PatchInfo(indPatch).f_R*PatchInfo(indPatch).dfdx_R.'-PatchInfo(indPatch).f_L*PatchInfo(indPatch).dfdx_L.';
        Ki=(BoundaryTerm-Ki);
    elseif(strcmp(PatchInfo(indPatch).PolyType,'Laguerre')) %-- UP TO THIS MOMENT, NUMERICAL
        Mi=conj(PatchInfo(indPatch).f.*Weight)*PatchInfo(indPatch).f.';
        Ki=conj(PatchInfo(indPatch).dfdx.*Weight)*PatchInfo(indPatch).dfdx.';
    end
    Mtrue=Mtrue+PatchInfo(indPatch).coeff_mn'*Mi*PatchInfo(indPatch).coeff_mn;
    M=M+PatchInfo(indPatch).coeff_mn'*(PatchInfo(indPatch).kFlux.*Mi)*PatchInfo(indPatch).coeff_mn;
    K=K+PatchInfo(indPatch).coeff_mn'*Ki*PatchInfo(indPatch).coeff_mn;
end

[EigVc,EigVl]=eig(-K+k0^2*M);
EigVl=diag(EigVl);

indEigVlreal=find(imag(EigVl)==0);
EigVl=EigVl(indEigVlreal);
EigVc=EigVc(:,indEigVlreal);
[EigVl,Indexes]=sort(EigVl/k0^2,'ascend');
neff=sqrt(EigVl(end+1-indMode))

EigVc=EigVc(:,Indexes(end+1-indMode));

for indPatch=1:length(PatchInfo)

    x=PatchInfo(indPatch).x;
    f=EigVc.'*PatchInfo(indPatch).coeff_mn.'*PatchInfo(indPatch).f;

    figure(1)
    axis on
    grid on
    hold on
    plot(x,f,'LineWidth',2)
    xlabel('z (mm)')
    
end
