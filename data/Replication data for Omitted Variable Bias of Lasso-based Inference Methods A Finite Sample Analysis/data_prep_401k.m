%% Preliminaries
% Data preparation for the study of the effect of 401(k) eligibility
% Authors: Kaspar Wuthrich and Ying Zhu
% DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied. 
% Questions/error reports: kwuthrich@ucsd.edu

clear;

cd '/Users/kasparwuthrich/Dropbox/research/Kaspar_Ying_research/submission/ReStat/final submission/replication_package_final'

rng('default');

%% Prepare data for analysis in R (Code and data are adapted from the original replication package for Belloni et al. (2017,ECMA))

for Inspec=[1 2]
    
    spec = Inspec;
    disp(spec);
    
    %% Import data
    data = dlmread('restatw.dat','\t',1,0);
    wgt = data(:,1);
    ira	= data(:,2);
    a401 = data(:,3);
    hval = data(:,4);
    hmort = data(:,5);
    hequity	= data(:,6);
    nifa = data(:,7);
    net_nifa = data(:,8);
    tfa = data(:,9);
    net_tfa	= data(:,10);
    tfa_he = data(:,11);
    tw = data(:,12);
    age = data(:,13);
    inc = data(:,14);
    fsize = data(:,15);
    educ = data(:,16);
    db = data(:,17);
    marr = data(:,18);
    male = data(:,19);
    twoearn = data(:,20);
    dum91 = data(:,21);
    e401 = data(:,22);
%     i1 = data(:,23);
    i2 = data(:,24);
    i3 = data(:,25);
    i4 = data(:,26);
    i5 = data(:,27);
    i6 = data(:,28);
    i7 = data(:,29);
    p401 = data(:,30);
    pira = data(:,31);
    hown = data(:,32);
    a1 = data(:,33);
    a2 = data(:,34);
    a3 = data(:,35);
    a4 = data(:,36);
    a5 = data(:,37);
    nohs = data(:,38);
    hs = data(:,39);
    smcol = data(:,40);
    col = data(:,41);
    a = data(:,42);
    icat = data(:,43);
    ecat = data(:,44);
    f1 = data(:,45);
    f2 = data(:,46);
    f3 = data(:,47);
    f4 = data(:,48);
    f5 = data(:,49);
    f6 = data(:,50);
    f = data(:,51);
    zhat = data(:,52);
    wntfa1 = data(:,53);
    wntfa2 = data(:,54);
    wntfa3 = data(:,55);
    wntfa4 = data(:,56);
    wntfa5 = data(:,57);
    wntfa6 = data(:,58);
    wntfa7 = data(:,59);
    wnet_nifa1 = data(:,60);
    wnet_nifa2 = data(:,61);
    wnet_nifa3 = data(:,62);
    wnet_nifa4 = data(:,63);
    wnet_nifa5 = data(:,64);
    wnet_nifa6 = data(:,65);
    wnet_nifa7 = data(:,66);
    wtw1 = data(:,67);
    wtw2 = data(:,68);
    wtw3 = data(:,69);
    wtw4 = data(:,70);
    wtw5 = data(:,71);
    wtw6 = data(:,72);
    wtw7 = data(:,73);
    net_n401 = data(:,74);
    wnet_n4014 = data(:,75);
    i1 = icat == 1;
    
    clear data;

    %% Save data for both specifications

    y = tw;
    d = e401;
    n = size(y,1);

    if spec == 1

        xMain = [i2,i3,i4,i5,i6,i7,a2,a3,a4,a5,fsize,hs,smcol,col,marr,twoearn,db,pira,hown]; % Benjamin, CH baseline

        pMain = size(xMain,2);
        xInt = zeros(n,pMain*pMain);
        for ii = 1:pMain
            for jj = 1:pMain
                xInt(:,pMain*(ii-1)+jj) = xMain(:,ii).*xMain(:,jj);
            end
        end
        x = [xMain xInt];
        keep = (std(x) > 1e-5);
        x = x(:,keep);
        [~,keepA] = noncollinear_Revision(x);

        Xmat = x(:,keepA);
        size(x,2) % p including collinear terms
        size(Xmat,2)  % p excluding collinear terms

        save('data_spec1');
        
    elseif spec == 2

        age = (age-25)/(64-25);
        inc = (inc+2652)/(242124+2652);
        fsize = fsize/13;
        educ = educ/18;
        idum = [i1,i2,i3,i4,i5,i6,i7];
        adum = [a1,a2,a3,a4,a5];
        edum = [nohs hs smcol col];
        odum = [marr twoearn db pira hown];
        iS = [idum inc (inc*ones(1,7)).*idum inc.^2 ((inc.^2)*ones(1,7)).*idum];
        aS = [age age.^2 age.^3];
        eS = [educ educ.^2];
        fS = [fsize fsize.^2];
        xMain = [iS aS eS fS odum];
        xO = [aS eS fS odum];
        nI = size(iS,2);
        nO = size(xO,2);
        xInt = zeros(n,nI*nO);
        for ii = 1:nI
            for jj = 1:nO
                xInt(:,nO*(ii-1)+jj) = iS(:,ii).*xO(:,jj);
            end
        end
        x = [xMain xInt];
        keep = (std(x) > 1e-5);
        x = x(:,keep);
        [~,keepA] = noncollinear_Revision(x);

        Xmat = x(:,keepA);
        size(x,2) % p including collinear terms
        size(Xmat,2)  % p excluding collinear terms

        save('data_spec2');
                
    end

end    