
for i = 1:40
    tst = strcat('s',int2str(i));
    imgPathb = 'C:/Users/Chetan There/Desktop/images/';
    imgPath = strcat(imgPathb,tst);
    imgPath
    imgType = '/*.pgm';
    images  = dir([imgPath imgType]);
    l = length(images) ;
    k = 0;
    for file = images'        
        fn = file.name;
        imgPathn = strcat(imgPath,'/',fn);
        imgPathn
        A = imread(imgPathn);
        %image(A)        
        A1 = reshape(A',1,10304);
        k = k + 1  ;     
        if k == 1
            %if tst == 's1'
            tmp = strcmp(tst,'s1');
            if tmp == 1
                xtrain = A1;
            else
                xtrain = [xtrain;A1];
            end
                        
        elseif k <= 5
            xtrain = [xtrain;A1];
            
        elseif k == 6
            
            tmp = strcmp(tst,'s1');
            if tmp == 1
                xtest = A1;
            else
                xtest = [xtest;A1];
            end                
            
        else
            xtest = [xtest;A1];
            
        end  
    end
    
end







%----------------------------------------------
%Step 1.1 KNN
%1. 5 fold cross validation




acounterl = [];
for j = 1:5
    kxtest  = [];
    kxtrain = [];
    
    k = j;
      
    for i = k:5:200
        A1 = xtest(i,:);       
        kxtest = [kxtest;A1];
              
        %add train data also
        A1     = xtrain(i,:);
        kxtest = [kxtest;A1];
        
        %forming train data for KNN
        ls = mod(i,5);
        l  = i - ls + 1;
        %changes for j = 5
        if j == 5
            l = l - 5;
        end
        
        %if tflag == 0
        for p = l:(l + 4)
            if p ~= i
                
                A1 = xtest(p,:);
                    
                kxtrain = [kxtrain;A1];
                    
                %add train data also
                A1 = xtrain(p,:);
                kxtrain = [kxtrain;A1];                
            end    
        end
       
    end    
    
    
           
    %Apply KNN here
    cc = 0;
    acounter = 0;
    
        for i = 1:2:80
            tl = i + 1;
            cc = cc + 1;
            for p = i:tl
                chkd = [];
                t = kxtest(p,:);
                %cal dist from all train ex
                for l = 1:320
                    tt = kxtrain(l,:);
                    d = 0;
                    
                    for  at = 1:10304
                        dat =  double(t(at)) - double(tt(at));
                        dat = dat * dat;
                          d = double(d) + double(dat);
                    end
                    d = abs(d);
                    chkd = [chkd d];
                    
                end
                %find min of 320
                schkd = size(chkd);
                
                [minv mini] = min(chkd);
               
                mdi = floor(mini / 8);
                mdf = (mini/8) - mdi;
                
                if mdf == 0
                    pc = mdi;
                else
                    pc = mdi + 1;
                end

                if pc == cc
                    acounter = acounter + 1;
                end

            end
        end
        acounterl = [acounterl acounter];
   
end

acc_sum = 0;
for i = 1:5
    acc_sum = acc_sum + acounterl(i);
end


acc_tot_knn = acc_sum * 100 / 400;

x = ['acc for KNN with 5 k fold = ', num2str(acc_tot_knn)];
disp(x)



%----------------------------------------------
%Step 1.2 KNN using PCA
%   5 fold cross validation



acounterl = [];
for j = 1:5
    kxtest  = [];
    kxtrain = [];
    
    k = j;
      
    for i = k:5:200
        A1 = xtest(i,:);       
        kxtest = [kxtest;A1];
              
        %add train data also
        A1     = xtrain(i,:);
        kxtest = [kxtest;A1];
        
        %forming train data for KNN
        ls = mod(i,5);
        l  = i - ls + 1;
        %changes for j = 5
        if j == 5
            l = l - 5;
        end
        
        %if tflag == 0
        for p = l:(l + 4)
            if p ~= i
                
                A1 = xtest(p,:);
                    
                kxtrain = [kxtrain;A1];
                    
                    %add train data also
                A1 = xtrain(p,:);
                kxtrain = [kxtrain;A1];                
            end    
        end
       
    end    
    
    
    
    %PCA
    pkxtest  = [];
    pkxtrain = [];
    trainpca = [];
    testpca  = [];
    
    
    
    pkxtrain = kxtrain';
    size(pkxtrain);
    [M N] = size(pkxtrain);
    mn = mean(pkxtrain,2);
    pkxtrain = double(pkxtrain) - repmat(mn,1,N);
   
    %covariance =  ( double(pkxtrain) * double(pkxtrain'));
    %[PC, V] = eig(covariance);
    
    covariance = cov(double(pkxtrain'));
    [PC, V] = eigs(covariance,80);  
       
    
    signals = PC' * double(pkxtrain);   
    ssig = size(signals);
    
   
    plot(V);
    

    pkxtest = kxtest';
    size(pkxtest);
    [M N] = size(pkxtest);
    mn = mean(pkxtest,2);
    pkxtest = double(pkxtest) - repmat(mn,1,N);
    
    
    
    covarianceT = cov(double(pkxtest'));
    [PCT, VT] = eigs(covarianceT,80); 
    testsignals = PCT' * double(pkxtest) ;
    ssigt = size(testsignals);


    trainpca = signals' ;
    testpca  = testsignals' ;

       
        
    %Apply KNN here using PCA
    cc = 0;
    acounter = 0;
    
    for i = 1:2:80
        tl = i + 1;
        cc = cc + 1;
        for p = i:tl
            chkd = [];
            t = testpca(p,:);
            %cal dist from all train ex
            for l = 1:320
                tt = trainpca(l,:);
                d = 0;

                for  at = 1:80
                    dat =  double(t(at)) - double(tt(at));
                    dat = dat * dat;
                      d = double(d) + double(dat);
                end
                d = abs(d);
                chkd = [chkd d];

            end
            %find min of 320
            schkd = size(chkd);
            chkd;
            [minv mini] = min(chkd);
            
            mdi = floor(mini / 8);
            mdf = (mini/8) - mdi;
            
            if mdf == 0
                pc = mdi;
            else
                pc = mdi + 1;
            end
            
            if pc == cc
                acounter = acounter + 1;
            end

        end
    end
    acounterl = [acounterl acounter];
   
end

acc_sum = 0;
for i = 1:5
    acc_sum = acc_sum + acounterl(i);
end


acc_tot_knn = acc_sum * 100 / 400

x = ['acc for KNN with 5 k fold and PCA  = ', num2str(acc_tot_knn)];
disp(x)





%Step 2:

   

%RESIZE IMAGE N TEST AGAIN

for i = 1:40
    tst = strcat('s',int2str(i));
    imgPathb = 'C:/Users/Chetan There/Desktop/images/';
    imgPath = strcat(imgPathb,tst);
    imgPath
    imgType = '/*.pgm';
    images  = dir([imgPath imgType]);
    l = length(images) ;
    k = 0;
    for file = images'        
        fn = file.name;
        imgPathn = strcat(imgPath,'/',fn);
        imgPathn
        A = imread(imgPathn);
        %image(A)        
        %A1 = reshape(A',1,10304);
        B = imresize(A, [56 46]);
        B1 = reshape(B',1,2576);
        k = k + 1  ;     
        if k == 1
            %if tst == 's1'
            tmp = strcmp(tst,'s1');
            if tmp == 1
                rxtrain = B1;
            else
                rxtrain = [rxtrain;B1];
            end
                        
        elseif k <= 5
            rxtrain = [rxtrain;B1];
            
        elseif k == 6
            
            tmp = strcmp(tst,'s1');
            if tmp == 1
                rxtest = B1;
            else
                rxtest = [rxtest;B1];
            end                
            
        else
            rxtest = [rxtest;B1];
            
        end  
    end
    
end

rxtests  = size(rxtest)
rxtrains = size(rxtrain)







%KNN ON RESIZE
acounterl = [];
for j = 1:5
    rkxtest  = [];
    rkxtrain = [];
    
    k = j;
      
    for i = k:5:200
        A1 = rxtest(i,:);
        
        rkxtest = [rkxtest;A1];
        
        
        %add train data also
        A1     = rxtrain(i,:);
        rkxtest = [rkxtest;A1];
        
        %forming train data for KNN
        ls = mod(i,5);
        l  = i - ls + 1;
        %changes for j = 5
        if j == 5
            l = l - 5;
        end
        
        for p = l:(l + 4)
            if p ~= i
                
                A1 = rxtest(p,:);
                    
                rkxtrain = [rkxtrain;A1];
                    
                A1 = rxtrain(p,:);
                rkxtrain = [rkxtrain;A1];                
            end    
        end
        %end
    end    
    
    
    %Apply KNN here ON RESIZE
    cc = 0;
    acounter = 0;
    %if tflag == 0
        for i = 1:2:80
            tl = i + 1;
            cc = cc + 1;
            for p = i:tl
                chkd = [];
                t = rkxtest(p,:);
                %cal dist from all train ex
                for l = 1:320
                    tt = rkxtrain(l,:);
                    d = 0;
                    %dat = int64(0);
                    for  at = 1:2576
                        dat =  double(t(at)) - double(tt(at));
                        dat = dat * dat;
                          d = double(d) + double(dat);
                    end
                    d = abs(d);
                    chkd = [chkd d];
                    
                end
                %find min of 320
                schkd = size(chkd);
                chkd;
                [minv mini] = min(chkd);
                i;
                mdi = floor(mini / 8);
                mdf = (mini/8) - mdi;
                
                if mdf == 0
                    pc = mdi;
                else
                    pc = mdi + 1;
                end

                if pc == cc
                    acounter = acounter + 1;
                end

            end
        end
        acounterl = [acounterl acounter];
    %end
end

acc_sum = 0;
for i = 1:5
    acc_sum = acc_sum + acounterl(i);
end

acc_tot_rknn = acc_sum * 100 / 400;

x = ['acc for resize and KNN with 5 k fold   = ', num2str(acc_tot_rknn)];
disp(x)

    


%KNN on resize using PCA

acounterl = [];
for j = 1:5
    rkxtest  = [];
    rkxtrain = [];
    
    k = j;
      
    for i = k:5:200
        A1 = rxtest(i,:);
        
        rkxtest = [rkxtest;A1];
        
        
        %add train data also
        A1     = rxtrain(i,:);
        rkxtest = [rkxtest;A1];
        
        %forming train data for KNN
        ls = mod(i,5);
        l  = i - ls + 1;
        %changes for j = 5
        if j == 5
            l = l - 5;
        end
        
        for p = l:(l + 4)
            if p ~= i
                
                A1 = rxtest(p,:);
                    
                rkxtrain = [rkxtrain;A1];
                    
                A1 = rxtrain(p,:);
                rkxtrain = [rkxtrain;A1];                
            end    
        end
        %end
    end    
    
    %PCA
    
    pkxtest  = [];
    pkxtrain = [];
    trainpca = [];
    testpca  = [];
       
    
    pkxtrain = rkxtrain';
    size(pkxtrain);
    [M N] = size(pkxtrain);
    mn = mean(pkxtrain,2);
    pkxtrain = double(pkxtrain) - repmat(mn,1,N);
   
    %covariance =  ( double(pkxtrain) * double(pkxtrain'));
    %[PC, V] = eig(covariance);
    
    covariance = cov(double(pkxtrain'));
    [PC, V] = eigs(covariance,80);    
    
           
    signals = PC' * double(pkxtrain);   
    ssig = size(signals);
    
   
       

    pkxtest = rkxtest';
    size(pkxtest);
    [M N] = size(pkxtest);
    mn = mean(pkxtest,2);
    pkxtest = double(pkxtest) - repmat(mn,1,N);
    
    
    
    covarianceT = cov(double(pkxtest'));
    [PCT, VT] = eigs(covarianceT,80); 
    testsignals = PCT' * double(pkxtest) ;
    ssigt = size(testsignals);


    trainpca = signals' ;
    testpca  = testsignals' ;
       
      
    
    
    
    %Apply KNN here ON RESIZE
    cc = 0;
    acounter = 0;
    %if tflag == 0
        for i = 1:2:80
            tl = i + 1;
            cc = cc + 1;
            for p = i:tl
                chkd = [];
                t = testpca(p,:);
                %cal dist from all train ex
                for l = 1:320
                    tt = trainpca(l,:);
                    d = 0;
                    %dat = int64(0);
                    for  at = 1:80
                        dat =  double(t(at)) - double(tt(at));
                        dat = dat * dat;
                          d = double(d) + double(dat);
                    end
                    d = abs(d);
                    chkd = [chkd d];
                    
                end
                %find min of 320
                schkd = size(chkd);
                chkd;
                [minv mini] = min(chkd);
                i;
                mdi = floor(mini / 8);
                mdf = (mini/8) - mdi;
                
                if mdf == 0
                    pc = mdi;
                else
                    pc = mdi + 1;
                end

                if pc == cc
                    acounter = acounter + 1;
                end

            end
        end
        acounterl = [acounterl acounter];
    %end
end

acc_sum = 0;
for i = 1:5
    acc_sum = acc_sum + acounterl(i);
end
acc_sum = 400 - acc_sum;
acc_tot_rknn = acc_sum * 100 / 400;

x = ['acc for resize and KNN with 5 k fold n PCA  = ', num2str(acc_tot_rknn)];
disp(x)

    




%Task 3:
%LDA


xtests  = size(xtest);
xtrains = size(xtrain);

acounterl = [];
for j = 1:5
    kxtest  = [];
    kxtrain = [];
    
    k = j;
      
    for i = k:5:200
        A1 = xtest(i,:);       
        kxtest = [kxtest;A1];
              
        %add train data also
        A1     = xtrain(i,:);
        kxtest = [kxtest;A1];
        
        %forming train data for KNN
        ls = mod(i,5);
        l  = i - ls + 1;
        %changes for j = 5
        if j == 5
            l = l - 5;
        end
        
        %if tflag == 0
        for p = l:(l + 4)
            if p ~= i
                
                A1 = xtest(p,:);                    
                kxtrain = [kxtrain;A1];
                    
                %add train data also
                A1 = xtrain(p,:);
                kxtrain = [kxtrain;A1];                
            end    
        end
       
    end    
    
    kxtrains  = size(kxtrain);
    kxtests   = size(kxtest);
    
   
    gmean = zeros(10304,1);
    sw = zeros(10304,10304);
    
    gmeantest = zeros(10304,1);
    swtest = zeros(10304,10304);
  
    for i = 0:39
        is = i * 8 + 1;
        ie = is + 7;

        tst = kxtrain(is:ie,:);            
        tst = tst';
        tst = double(tst);
        
        mn = mean(tst,2);        
     
        gmean = gmean + mn;        
        
        cv1 = tst - repmat(mn,1,8);
        cv = cv1 * cv1';        
        sw = sw + cv;
        %size(cv)             
        
        
        %FOR TEST DATA
        istest = i * 2 + 1;
        ietest = istest + 1;
        test = kxtest(istest:ietest,:);            
        test = test';
        test = double(test);
        mntest = mean(test,2); 
        gmeantest = gmeantest + mntest;   
        cv1test = test - repmat(mntest,1,2);
        cvtest = cv1test * cv1test';        
        swtest = swtest + cvtest;
        
    end
    
       
    sgmean = size(gmean);
    ssw = size(sw)    ;
    
    sgmeantest = size(gmeantest);
    sswtest = size(swtest)   ;
    sb = zeros(10304,10304);
    sbtest = zeros(10304,10304);
    
    for i = 0:39
        is = i * 8 + 1;
        ie = is + 7;

        in = i + 1;
        tst = kxtrain(is:ie,:);        
        tst = tst';
        tst = double(tst);
        %size(tst)
        mn = mean(tst,2);        
        mtst = (mn - gmean)* (mn - gmean)';        
        sb = sb + mtst;    
        
                   
        
        
        %FOR TEST DATA
        istest = i * 2 + 1;
        ietest = istest + 1;
        test = kxtest(istest:ietest,:);            
        test = test';
        test = double(test);
        mntest = mean(test,2); 
        mtest = (mntest - gmeantest)* (mntest - gmeantest)'; 
        sbtest = sbtest + mtest; 
        
    end
    
    
     ssb = size(sb);
     
     
     invsw = inv(sw);
     invswsb = invsw * sb;
     
     sswf = size(invswsb);
     
     [V, D] = eigs(invswsb,39);
     
     
     
     %final train data set     
     ldatrain = double(kxtrain) * V;     
     sldatrain = size(ldatrain);
     
    
     
     %FOR TEST DATA
     
     ssbtest = size(sbtest);
     
     invswtest = inv(swtest);
     invswsbtest = invswtest * sbtest;
     
     sinvswsbtest = size(invswsbtest);
     
     [VTEST, DTEST] = eigs(invswsbtest,39);
     
     
     
     %final train data set     
     ldatest = double(kxtest) * VTEST;     
     sldatest = size(ldatest);
     
     
        
    %Apply KNN here
    cc = 0;
    acounter = 0;
    
        for i = 1:2:80
            tl = i + 1;
            cc = cc + 1;
            for p = i:tl
                chkd = [];
                t = ldatest(p,:);
                %cal dist from all train ex
                for l = 1:320
                    tt = ldatrain(l,:);
                    d = 0;
                    
                    for  at = 1:39
                        dat =  double(t(at)) - double(tt(at));
                        dat = dat * dat;
                          d = double(d) + double(dat);
                    end
                    d = abs(d);
                    chkd = [chkd d];
                    
                end
                %find min of 320
                schkd = size(chkd);
                chkd;
                [minv mini] = min(chkd);
                i;
                mdi = floor(mini / 8);
                mdf = (mini/8) - mdi;
                mini;
                cc;
                if mdf == 0
                    pc = mdi;
                else
                    pc = mdi + 1;
                end

                if pc == cc
                    acounter = acounter + 1;
                end

            end
        end
        acounterl = [acounterl acounter];
   
end
acounterl;
acc_sum = 0;
for i = 1:5
    acc_sum = acc_sum + acounterl(i);
end
acc_sum = 400 - acc_sum;

acc_tot_knn = acc_sum * 100 / 400

x = ['acc for KNN with 5 k fold and LDA  = ', num2str(acc_tot_knn)];
disp(x)







%Task 4: 
%PCA and LDA and do KNN


acounterl = [];
for j = 1:5
    kxtest  = [];
    kxtrain = [];
    
    k = j;
      
    for i = k:5:200
        A1 = xtest(i,:);       
        kxtest = [kxtest;A1];
              
        %add train data also
        A1     = xtrain(i,:);
        kxtest = [kxtest;A1];
        
        %forming train data for KNN
        ls = mod(i,5);
        l  = i - ls + 1;
        %changes for j = 5
        if j == 5
            l = l - 5;
        end
        
        %if tflag == 0
        for p = l:(l + 4)
            if p ~= i
                
                A1 = xtest(p,:);
                    
                kxtrain = [kxtrain;A1];
                    
                    %add train data also
                A1 = xtrain(p,:);
                kxtrain = [kxtrain;A1];                
            end    
        end
       
    end    
    
    kxtrains  = size(kxtrain);
    kxtests   = size(kxtest);
    
    %PCA
    pkxtest  = [];
    pkxtrain = [];
    trainpca = [];
    testpca  = [];
    
    
    
    pkxtrain = kxtrain';
    size(pkxtrain);
    [M N] = size(pkxtrain)
    mn = mean(pkxtrain,2);
    pkxtrain = double(pkxtrain) - repmat(mn,1,N);
    
    
    covariance = cov(double(pkxtrain'));
    [PC, V] = eigs(covariance,80);  
    
    
    
    signals = PC' * double(pkxtrain);
    
    %signals = pkxtrain * PC';
    ssig = size(signals);
    
    

    pkxtest = kxtest';
    size(pkxtest);
    [M N] = size(pkxtest);
    mn = mean(pkxtest,2);
    pkxtest = double(pkxtest) - repmat(mn,1,N);
    
    
    
    covarianceT = cov(double(pkxtest'));
    [PCT, VT] = eigs(covarianceT,80); 
    testsignals = PCT' * double(pkxtest) ;
    ssigt = size(testsignals);


    trainpca = signals' ;
    testpca  = testsignals' ;

    spca = size(trainpca);
    stpca = size(testpca);
     
    %Apply LDA on PCA now
    
    gmean = zeros(80,1);
    sw = zeros(80,80);
    
    gmeantest = zeros(80,1);
    swtest = zeros(80,80);
  
    for i = 0:39
        is = i * 8 + 1;
        ie = is + 7;

        tst = trainpca(is:ie,:);            
        tst = tst';
        tst = double(tst);
        %size(tst)
        mn = mean(tst,2);        
        %size(mn)
        gmean = gmean + mn;        
        %tst = tst - repmat(mn,1,8);        
        %cv = cov(tst');
        cv1 = tst - repmat(mn,1,8);
        cv = cv1 * cv1';  
        size(cv)  ;
        sw = sw + cv;
            
        
        
        %FOR TEST DATA
        istest = i * 2 + 1;
        ietest = istest + 1;
        test = testpca(istest:ietest,:);            
        test = test';
        test = double(test);
        mntest = mean(test,2); 
        gmeantest = gmeantest + mntest;   
        cv1test = test - repmat(mntest,1,2);
        cvtest = cv1test * cv1test';        
        swtest = swtest + cvtest;
        
    end
    
       
    sgmean = size(gmean);
    ssw = size(sw)  ;  
    
    sgmeantest = size(gmeantest);
    sswtest = size(swtest)   ;
    sb = zeros(80,80);
    sbtest = zeros(80,80);
    
    for i = 0:39
        is = i * 8 + 1;
        ie = is + 7;

        in = i + 1;
        tst = trainpca(is:ie,:);        
        tst = tst';
        tst = double(tst);
        %size(tst)
        mn = mean(tst,2);        
        mtst = (mn - gmean)* (mn - gmean)';        
        %size(mtst)        
        %size(mtst')
        sb = sb + mtst;        
        %cv = cov(mtst');
        %sb = sb + cv;
        %size(cv)             
        
        
        %FOR TEST DATA
        istest = i * 2 + 1;
        ietest = istest + 1;
        test = testpca(istest:ietest,:);            
        test = test';
        test = double(test);
        mntest = mean(test,2); 
        mtest = (mntest - gmeantest)* (mntest - gmeantest)'; 
        sbtest = sbtest + mtest; 
        
    end
    
    
     ssb = size(sb);
     
     
     invsw = inv(sw);
     invswsb = invsw * sb;
     
     sswf = size(invswsb);
     
     [V, D] = eigs(invswsb,39);
     
  ;
     
     %final train data set     
     ldatrain = double(trainpca) * V;     
     sldatrain = size(ldatrain);
     
    
     
     %FOR TEST DATA
     
     ssbtest = size(sbtest);
     
     invswtest = inv(swtest);
     invswsbtest = invswtest * sbtest;
     
     sinvswsbtest = size(invswsbtest)
     
     [VTEST, DTEST] = eigs(invswsbtest,39);
     
     
     sVTEST = size(VTEST);
     
     %final train data set     
     ldatest = double(testpca) * VTEST;     
     sldatest = size(ldatest);    
    
    
     
    
        
    %Apply KNN here using PCA
    cc = 0;
    acounter = 0;
    
    for i = 1:2:80
        tl = i + 1;
        cc = cc + 1;
        for p = i:tl
            chkd = [];
            t = ldatest(p,:);
            %cal dist from all train ex
            for l = 1:320
                tt = ldatrain(l,:);
                d = 0;

                for  at = 1:39
                    dat =  double(t(at)) - double(tt(at));
                    dat = dat * dat;
                      d = double(d) + double(dat);
                end
                d = abs(d);
                chkd = [chkd d];

            end
            %find min of 320
            schkd = size(chkd);
            chkd;
            [minv mini] = min(chkd);
            i;
            mdi = floor(mini / 8);
            mdf = (mini/8) - mdi;
            
            if mdf == 0
                pc = mdi;
            else
                pc = mdi + 1;
            end
            
            if pc == cc
                acounter = acounter + 1;
            end

        end
    end
    acounterl = [acounterl acounter];
   
end
acounterl;
acc_sum = 0;
for i = 1:5
    acc_sum = acc_sum + acounterl(i);
end
acc_sum = 400 - acc_sum;


acc_tot_knn = acc_sum * 100 / 400

x = ['acc for KNN with 5 k fold PCA n LDA  = ', num2str(acc_tot_knn)];
disp(x)






%Task 5:
%SVM using 5 fold cross validation



%SVM USING K FOLD
xtests  = size(xtest);
xtrains = size(xtrain);

acounterl = [];
for j = 1:5
    kxtest  = [];
    kxtrain = [];
    
    k = j;
      
    for i = k:5:200
        A1 = xtest(i,:);
        
        kxtest = [kxtest;A1];
        
        
        %add train data also
        A1     = xtrain(i,:);
        kxtest = [kxtest;A1];
        
        %forming train data for KNN
        ls = mod(i,5);
        l  = i - ls + 1;
        %changes for j = 5
        if j == 5
            l = l - 5;
        end
        
        for p = l:(l + 4)
            if p ~= i
                
                A1 = xtest(p,:);
                    
                kxtrain = [kxtrain;A1];
                    
                A1 = xtrain(p,:);
                kxtrain = [kxtrain;A1];                
            end    
        end
       
    end    
    
    kxtrains  = size(kxtrain);
    kxtests   = size(kxtest);
    
    
    %FORM HYPERPLANES HERE
    hypw0 = [];
    
    %
    hypw = []
    Z = []
    
    for i1 =  1:8:320
        j1 = i1 - 1;
        z1 = -ones(j1,1);
        z2 = ones(8,1);
        l1 = 320 - j1 - 8;
        z3 = -ones(l1,1);
        Z = [z1;z2;z3];

        H = (double(kxtrain) * double(kxtrain')).*(Z * Z') ;      
        f = -ones(320,1);    
        A = -eye(320);    
        a = zeros(320,1);    
        B = [[Z'];[zeros(319,320)]];   
        b = zeros(320,1) ;   
        lb= zeros(320,1);
        ub = ones(320,1) * 100;

        ax = quadprog(H+eye(320)*0.001,f,A,a,B,b);
        sax = size(ax);

        w1 = (ax.*Z)';    
        w2 = double(w1) * double(kxtrain);
        w = w2';



        w01 = 1 ;   
        x11 = kxtrain(i1,:);
        w02 = double(w') * double(x11');
        w0 = w01 - w02;
        %add slack of unit 1
        %w0 = w0 + 1;

        %FORM ARRAY OF HYPERPLANES
        wa = w';
        if i1 == 1
            hypw = wa;
        else
            hypw = [hypw;wa];
        end

        hypw0 = [hypw0 w0];


    end

    shyperplane = size(hypw);
    
    
    %Apply svm here
    cc = 0;
    acounter = 0;
    acc_counter = 0;
    wacc_counter = 0;
    %if tflag == 0
        for i = 1:2:80
            tl = i + 1;
            cc = cc + 1;
            for p = i:tl
                chkval = [];
                t = kxtest(p,:);
                %cal dist from all train ex
                for u = 1:40            
                    wtt  = hypw(u,:);
                    wt = wtt';
                    w0t = hypw0(u);
                    swt = size(wt);
                    sw0t = size(w0t);
                    chk = dot(double(wt),double(t))+ w0t; 
                    chkval = [chkval chk];
                end   
                %find max of 40
                %chkval
                schkval = size(chkval);
                [maxv maxi] = max(chkval);

            
                if maxi == cc
                    acc_counter = acc_counter + 1;
                elseif maxv == chkval(cc)
                    acc_counter = acc_counter + 1;           
                else
                    wacc_counter = wacc_counter + 1;
                end         
              
            end
        end
        
        
        acounterl = [acounterl acc_counter];
    %end
end
acounterl;
acc_sum = 0;
for i = 1:5
    acc_sum = acc_sum + acounterl(i);
end

acc_sum;
acc_tot_knn = acc_sum * 100 / 400

x = ['acc for SVM with 5 k fold   = ', num2str(acc_tot_knn)];
disp(x)





%Task 6:
%SVM using PCA and 5 fold cross validation


xtests  = size(xtest);
xtrains = size(xtrain);

acounterl = [];
for j = 1:5
    kxtest  = [];
    kxtrain = [];
    
    k = j;
      
    for i = k:5:200
        A1 = xtest(i,:);
        
        kxtest = [kxtest;A1];
       
        
        %add train data also
        A1     = xtrain(i,:);
        kxtest = [kxtest;A1];
        
        %forming train data for KNN
        ls = mod(i,5);
        l  = i - ls + 1;
        %changes for j = 5
        if j == 5
            l = l - 5;
        end
        
        for p = l:(l + 4)
            if p ~= i
                
                A1 = xtest(p,:);
                    
                kxtrain = [kxtrain;A1];
                    
                A1 = xtrain(p,:);
                kxtrain = [kxtrain;A1];                
            end    
        end
        %end
    end    
    
    kxtrains  = size(kxtrain);
    kxtests   = size(kxtest);
    
    
    
    %PCA
    pkxtest  = [];
    pkxtrain = [];
    trainpca = [];
    testpca  = [];
    
    
    
    pkxtrain = kxtrain';
    size(pkxtrain);
    [M N] = size(pkxtrain);
    mn = mean(pkxtrain,2);
    pkxtrain = double(pkxtrain) - repmat(mn,1,N);
    
    covariance = cov(double(pkxtrain'));
    [PC, V] = eigs(covariance,80);  
    %[PC, V] = eig(covariance);
    
         
      
    pcs = size(PC);
    pts = size(pkxtrain);
    signals = PC' * double(pkxtrain);
    
    %signals = pkxtrain * PC';
    ssig = size(signals);
   

    pkxtest = kxtest';
    size(pkxtest);
    [M N] = size(pkxtest);
    mn = mean(pkxtest,2);
    pkxtest = double(pkxtest) - repmat(mn,1,N);
    
        
    covarianceT = cov(double(pkxtest'));
    [PCT, VT] = eigs(covarianceT,80); 
    testsignals = PCT' * double(pkxtest) ;
    ssigt = size(testsignals);


    trainpca = signals' ;
    testpca  = testsignals' ;

    spca = size(trainpca);
    stpca = size(testpca);
    
    
    
    
    
    
    %FORM HYPERPLANES HERE
    hypw0 = [];
    
    %
    hypw = []
    Z = []
    
    for i1 =  1:8:320
        j1 = i1 - 1;
        z1 = -ones(j1,1);
        z2 = ones(8,1);
        l1 = 320 - j1 - 8;
        z3 = -ones(l1,1);
        Z = [z1;z2;z3];

        H = (double(trainpca) * double(trainpca')).*(Z * Z') ;      
        f = -ones(320,1);    
        A = -eye(320);    
        a = zeros(320,1);    
        B = [[Z'];[zeros(319,320)]];   
        b = zeros(320,1) ;   
        lb= zeros(320,1);
        ub = ones(320,1) * 100;

        ax = quadprog(H+eye(320)*0.001,f,A,a,B,b);
        sax = size(ax);

        w1 = (ax.*Z)';    
        w2 = double(w1) * double(trainpca);
        w = w2';



        w01 = 1 ;   
        x11 = trainpca(i1,:);
        w02 = double(w') * double(x11');
        w0 = w01 - w02;
        %add slack of unit 1
        %w0 = w0 + 1;

        %FORM ARRAY OF HYPERPLANES
        wa = w';
        if i1 == 1
            hypw = wa;
        else
            hypw = [hypw;wa];
        end

        hypw0 = [hypw0 w0];


    end

    shyperplane = size(hypw);
    
    
    %Apply svm here
    cc = 0;
    acounter = 0;
    acc_counter = 0;
    wacc_counter = 0;
    %if tflag == 0
        for i = 1:2:80
            tl = i + 1;
            cc = cc + 1;
            for p = i:tl
                chkval = [];
                t = testpca(p,:);
                %cal dist from all train ex
                for u = 1:40            
                    wtt  = hypw(u,:);
                    wt = wtt';
                    w0t = hypw0(u);
                    swt = size(wt);
                    sw0t = size(w0t);
                    chk = dot(double(wt),double(t))+ w0t; 
                    chkval = [chkval chk];
                end   
                %find max of 40
                %chkval
                schkval = size(chkval);
                [maxv maxi] = max(chkval);

            
                if maxi == cc
                    acc_counter = acc_counter + 1;
                elseif maxv == chkval(cc)
                    acc_counter = acc_counter + 1;           
                else
                    wacc_counter = wacc_counter + 1;
                end         
              
            end
        end
        
        
        acounterl = [acounterl acc_counter];
    %end
end
acounterl;
acc_sum = 0;
for i = 1:5
    acc_sum = acc_sum + acounterl(i);
end
acc_sum = 400 - acc_sum;

acc_tot_knn = acc_sum * 100 / 400;

x = ['acc for SVM with 5 k fold n PCA  = ', num2str(acc_tot_knn)];
disp(x)












%{

%PCA Test

kxtrains  = size(kxtrain)
kxtests   = size(kxtest)


pkxtrain = kxtrain';

[M N] = size(pkxtrain)

mnv = [];
for i1=1:M
    mn = mean(pkxtrain(i1,:));
    mnv = [mnv;mn];
end
%mn = mean(pkxtrain,2);
%pkxtrain = double(pkxtrain) - repmat(mn,1,N);
pkxtrain = double(pkxtrain) - repmat(mnv,1,N);
%}



%{
mnv = [];
for i1=1:M
    mn = mean(pkxtrain(i1,:));
    mnv = [mnv;mn];
end
mnvs = size(mnv)
pkxtrain2 = double(pkxtrain) - repmat(mnv,1,N);
covariance =  ( pkxtrain2 * pkxtrain2');
pkt = pkxtrain2(1,:);
%}

%{
covariance =  ( double(pkxtrain) * double(pkxtrain'));
[PC, V] = eig(covariance);
V = diag(V);
[junk, rindices] = sort(-1*V);


rindicesn = rindices((1:50),:)



%PC2 = PC(:,rindicesn) ;
%PC = PC(rindicesn,:) ;

PC = PC(:,rindicesn) ;
signals = PC' * pkxtrain;

%}

%{




pkxtest = kxtest';
size(pkxtest);
[M N] = size(pkxtest)



%mn = mean(pkxtrain,2);
%pkxtrain = double(pkxtrain) - repmat(mn,1,N);

mnv = [];
for i1=1:M
    mn = mean(pkxtest(i1,:));
    mnv = [mnv;mn];
end
mnvs = size(mnv)
pkxtest2 = double(pkxtest) - repmat(mnv,1,N);
covarianceT =  ( pkxtest2 * pkxtest2');
pkt = pkxtest2(1,:);

%covariance =  ( double(pkxtrain) * double(pkxtrain'));
[PCT, VT] = eig(covarianceT);
VT = diag(VT);
[junkT, rindicesT] = sort(-1*VT);

cvt = covariance(1,:) ;



V1T = VT(rindicesT);
%
rindicesnT = rindicesT((1:50),:);
%
V2T = VT(rindicesnT);

PC1T = PCT(:,rindicesT) ;
PC2T = PCT(:,rindicesnT) ;

size(VT)
size(PCT)
pct = PCT(1,:) ;
vt = VT;
signalsT = PCT2' * pkxtest2;

st = signals(1,:)







covariance = cov(double(pkxtrain'));

cv = covariance(1,:)
[PC, V] = eigs(covariance,80); 

PC
V

 pcs = size(PC);
 pts = size(pkxtrain);
 signals = PC' * double(pkxtrain);   
 
 st = signals(1,:)
 ssig = size(signals);




%signals = pkxtrain * PC;
ssig = size(signals)

%}



%{
pkxtest = kxtest';
size(pkxtest);
[M N] = size(pkxtest)
%mn = mean(pkxtest,2);
mnv = [];
for i1=1:M
    mn = mean(pkxtest(i1,:));
    mnv = [mnv;mn];
end
%pkxtest = double(pkxtest) - repmat(mn,1,N);
pkxtest = double(pkxtest) - repmat(mnv,1,N);
covarianceT = ( pkxtest * pkxtest');
%covarianceT = 1 / (N-1) * ( pkxtest' * pkxtest);


size(covarianceT)
[PCT, VT] = eig(covarianceT);
VT = diag(VT);
[junk, rindices] = sort(-1*VT);
VT = VT(rindices);
%PCT = PCT(:,rindices) ;
rindicesn = rindices((1:50),:);
%
VT = VT(rindicesn);

PCT = PCT(:,rindicesn) ;
%PCT = PCT(rindicesn,:) ;
size(VT)
size(PCT)
testsignals = PCT' * pkxtest;
%testsignals = PCT * pkxtest ;
ssigt = size(testsignals)


trainpca = signals' ;
testpca  = testsignals' ;

spca = size(trainpca)
stpca = size(testpca)
   
acounterl = [];
%PCA Test
%Apply KNN here
    cc = 0;
    acounter = 0;
    %if tflag == 0
        for i = 1:2:80
            tl = i + 1;
            cc = cc + 1;
            for p = i:tl
                chkd = [];
                t = testpca(p,:);
                %cal dist from all train ex
                for l = 1:320
                    tt = trainpca(l,:);
                    d = 0;
                    %dat = int64(0);
                    for  at = 1:50
                        dat =  double(t(at)) - double(tt(at));
                        dat = dat * dat;
                          d = double(d) + double(dat);
                    end
                    d = abs(d);
                    chkd = [chkd d];
                    %x = input('hi')
                end
                %find min of 320
                schkd = size(chkd);
                chkd;
                [minv mini] = min(chkd);
                i;
                mdi = floor(mini / 8);
                mdf = (mini/8) - mdi;
                mini;
                cc;
                if mdf == 0
                    pc = mdi;
                else
                    pc = mdi + 1;
                end

                if pc == cc
                    acounter = acounter + 1;
                end

            end
        end
        acounterl = [acounterl acounter]
    %end

    acounterl

 
 %}   



    

%RESIZE IMAGE N TEST AGAIN

for i = 1:40
    tst = strcat('s',int2str(i));
    imgPathb = 'C:/Users/Chetan There/Desktop/images/';
    imgPath = strcat(imgPathb,tst);
    imgPath
    imgType = '/*.pgm';
    images  = dir([imgPath imgType]);
    l = length(images) ;
    k = 0;
    for file = images'        
        fn = file.name;
        imgPathn = strcat(imgPath,'/',fn);
        imgPathn
        A = imread(imgPathn);
        %image(A)        
        %A1 = reshape(A',1,10304);
        B = imresize(A, [56 46]);
        B1 = reshape(B',1,2576);
        k = k + 1  ;     
        if k == 1
            %if tst == 's1'
            tmp = strcmp(tst,'s1');
            if tmp == 1
                rxtrain = B1;
            else
                rxtrain = [rxtrain;B1];
            end
                        
        elseif k <= 5
            rxtrain = [rxtrain;B1];
            
        elseif k == 6
            
            tmp = strcmp(tst,'s1');
            if tmp == 1
                rxtest = B1;
            else
                rxtest = [rxtest;B1];
            end                
            
        else
            rxtest = [rxtest;B1];
            
        end  
    end
    
end

rxtests  = size(rxtest)
rxtrains = size(rxtrain)



    


