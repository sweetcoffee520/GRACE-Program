function [DeltaGC,DeltaGS,Data_number,Data_time]=readdata(Data_type,Data_adress,Lmax,C20_address,C20_Flag)

%%Set the data provider and the data Address
DATA_Provider=Data_type;
Address=Data_adress;
%%------------------------------------------------------
%%You may check the max degree of the SH coefficients Lmax.
if isequal('CSR',DATA_Provider)
    Data_time=GrabTheDate_CSR(Address);
    [DeltaGC,DeltaGS,k]=read_GRACE_CSR(Address,Lmax,C20_address,C20_Flag);
%     Data_time=GrabTheDate_CSR(Address);
    disp('Data CSR is ready!')
else if isequal('GFZ',DATA_Provider)
        [DeltaGC,DeltaGS,k]=read_GRACE_GFZ(Address,Lmax);
        Data_time=GrabTheDate(Address);
        disp('Data GFZ is ready!')
    else if isequal('JPL',DATA_Provider)
            [DeltaGC,DeltaGS,k]=read_GRACE_JPL(Address,Lmax,C20_address,C20_Flag); 
            Data_time=GrabTheDate(Address);
            disp('Data JPL is ready!')
        else if isequal('GRGS',DATA_Provider)
                [DeltaGC,DeltaGS,k]=read_GRACE_GRGS(Address);
                Data_time=GrabTheDate_GRGS(Address);
                disp('Data GRGS is ready!')
               else if isequal('GLDAS',DATA_Provider)
                       [DeltaGC,DeltaGS,k]=read_GLDAS(Address,Lmax); 
                       Data_time=GrabTheDate_GLDAS(Address);
                       disp('Data GLDAS is ready!')
                   end
            end
        end
    end
end
Data_number=k;