function Iw = fog(I)

    I=double(I)/255;

    [row,col,z] = size(I);
    landline = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Iw = I;
    A = 0.8;
    m = floor(row/2);
    n = floor(col/2);
    for beta=0:0.01:0.1
        for i=1:z
            for j=landline+1:row
                for l=1:col
                    d(j,l) = 1/((j-landline)^.05 + 0.0001);
                    d2(j,l) = d(j,l)*8;
                    d(j,l) = -0.04*sqrt((j-m).^2+(l-n).^2) + 17;
                    td(j,l) = exp(-beta*d(j,l));
                    Iw(j,l,i) = I(j,l,i)*td(j,l) + A*(1-td(j,l));
                end
            end
        end
        str1 = ['D:\Disk_G\PROJECTS\Registration_baseline\testdata\',num2str(beta),'_fog.png'];
        imwrite(Iw,str1);
    end
end


