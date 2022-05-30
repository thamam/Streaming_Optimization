  %Try to use simulation samples from previous runs
    if UseSavedNhppSimdData
        try
            load(params.fileNhppName);
        catch ME
            if( strcmp(ME.identifier,'MATLAB:load:couldNotReadFile') )
                ME.message
                fprintf('Simulating new process \n');
                [alphaSyn] = syntContIntNHPP(params); %generate coeffs for Lambda(t)
                [ tK ] = nhppSynt(@(t)lambda_t(t,alphaSyn,params.phi_n,params.nInd) ,params); %send function handle of lambda as function of t
                dateOfFile=date;
                save(params.fileNhppName,'tK','alphaSyn','dateOfFile','params');
            else
                load(params.fileNhppName);
            end
        end
    else % Required new simulation
        fprintf('Simulating new process data \n');
        [alphaSyn] = syntContIntNHPP(params); %generate coeffs for Lambda(t)
        [ tK ] = nhppSynt(@(t)lambda_t(t,alphaSyn,params.phi_n,params.nInd) ,params); %send function handle of lambda as function of t
        if(OverWrite_NhppSimData)
            dateOfFile=date;
            save(params.fileNhppName,'tK','alphaSyn','dateOfFile');
        end
    end