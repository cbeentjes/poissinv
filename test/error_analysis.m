%%
% Plot the validation data derived from running poissinv_check.c and
% poissinv_check.cu for both single and double precision.
%

clear all
close all

check = [0,0,0,0];

%% Serial data
data_dir = 'Serial/';
file_name = 'poissinv_check.txt';

fid = fopen([data_dir,file_name]);
if fid ~= -1

    line = fgetl(fid);
    for kx = 1:3   
        % Skip empty lines
        for ix = 1:2
            line = fgetl(fid);
        end
        if line == -1
            continue
        end
        % Detect type of test
        switch line(1:6)
            case 'scalar'
                px = 1;
            case 'vector'
                px = 2;
            case '******'
                px = 3;
                for ix = 1:5
                    line = fgetl(fid);
                end
        end
        % Skip empty lines    
        for ix = 1:2
            line = fgetl(fid);
        end

        % Extract data
        line = fgetl(fid);
        l = 1;
        while length(line) > 10
           data = sscanf(line, '%f');
           lam(l,px) = data(1);
           E1_32(l,px) = data(2);
           E1_64(l,px) = data(3);

           line = fgetl(fid);
           l = l + 1;
        end  
        check(px) = 1;
    end
    fclose(fid);    
end

%% CUDA data
data_dir = 'CUDA/';
file_name = 'poissinv_check.txt';

fid = fopen([data_dir,file_name]);
if fid ~= -1
    line = fgetl(fid);
    px = 4;
    
    for ix = 1:3
        line = fgetl(fid);
    end
    line = fgetl(fid);
    l = 1;
    while length(line) > 10

       data = sscanf(line, '%f');
       lam(l,px) = data(1);
       E1_32(l,px) = data(2);
       E1_64(l,px) = data(3);

       line = fgetl(fid);
       l = l + 1;
       check(px) = 1;
    end  
    fclose(fid);        
end

% Fill in missing data
lam(:,check==0) = nan;
E1_32(:,check==0) = nan;
E1_64(:,check==0) = nan;

%% Plotting
close all

lg  = {'poissinv (CPU)' ,'poissinv\_v (CPU)' ,'poisscinv (CPU)', 'poissinv (GPU)'};
lgf = {'poissinvf (CPU)','poissinvf\_v (CPU)','poisscinvf (CPU)','poissinvf (GPU)'};

IX = 1:4;
IX = IX(check~=0);


figure(1)
for px = IX
    hold on
    plot(lam(:,px), E1_64(:,px),'.')
end
hold on
x = logspace(-4,7,1e3);
% Lam > 13.5
ix = x > 13.5;
plot(x(ix),3e-16*sqrt(x(ix)),'k')
% Lam < 13.5
ix = x <= 13.5;
plot(x(ix),1.6e-15*sqrt(max(1e-2,x(ix))),'k')

if intersect(IX,[2,4])
    % Vector version
    ix = x <= 4.5;
    plot(x(ix),1.2e-15*sqrt(max(1e-2,x(ix))),'b')
    ix = x > 4.5;
    plot(x(ix),2e-16*sqrt(max(160,x(ix))),'b')
end

xlabel('$\lambda$','interpreter','latex')
ylabel('$L_1$ error','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'YScale','log')
set(gca,'XScale','log')
legend(lg(IX),'location','best')
title('Double precision')

figure(2)
for px = IX
    hold on
    plot(lam(:,px), E1_32(:,px),'.')
end
hold on
x = logspace(-4,7,1e3);

% Lam > 11.5
ix = x > 11.5;
plot(x(ix),1.8e-7*sqrt(x(ix)),'k')
% Lam < 11.5
ix = x <= 11.5;
plot(x(ix),6.8e-7*sqrt(max(1e-2,x(ix))),'k')

if any(IX==2)
    % Vector CPU version
    ix = x > 6;
    plot(x(ix),1.8e-7*sqrt(max(69,x(ix))),'b')
    
    ix = x <= 6;
    plot(x(ix),6.2e-7*sqrt(max(1e-2,x(ix))),'b')
end

if any(IX==4)
    % Vector GPU version
    ix = x > 7.5;
    plot(x(ix),1.8e-7*sqrt(max(89,x(ix))),'r')
    
    ix = x <= 7.5;
    plot(x(ix),6.2e-7*sqrt(max(1e-2,x(ix))),'r')
end

xlabel('$\lambda$','interpreter','latex')
ylabel('$L_1$ error','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'YScale','log')
set(gca,'XScale','log')
legend(lg(IX),'location','best')
title('Single precision')
