function objFuns = ObjectiveFunctionStructure()

    %Forrester Function
    forr = struct();
    forr.Name = 'Forrester';
    forr.Func = @Forrester;
    forr.MinVals = [0];
    forr.MaxVals = [1];
    
    %Sphere Function
    sphere = struct();
    sphere.Name = 'Sphere';
    sphere.Func = @Sphere;
    sphere.MinVals = [-5.12, -5.12];
    sphere.MaxVals = [5.12, 5.12];
    
    %Booth Function
    booth = struct();
    booth.Name = 'Booth';
    booth.Func = @Booth;
    booth.MinVals = [-10, -10];
    booth.MaxVals = [10, 10];
    
    %Rosenbrock Function
    rosenbrock = struct();
    rosenbrock.Name = 'Rosenbrock';
    rosenbrock.Func = @Rosenbrock;
    rosenbrock.MinVals = [-2.408, -2.408];
    rosenbrock.MaxVals = [2.408, 2.408];
    
    %Beale Function
    beale = struct();
    beale.Name = 'Beale';
    beale.Func = @Beale;
    beale.MinVals = [-4.5, -4.5];
    beale.MaxVals = [4.5, 4.5];
    
    %GoldsteinPrice Function
    gp = struct();
    gp.Name = 'GoldsteinPrice';
    gp.Func = @GoldsteinPrice;
    gp.MinVals = [-2,-2];
    gp.MaxVals = [2, 2];
    
    %Zakharov Function
    zak = struct();
    zak.Name = 'Zakharov';
    zak.Func = @Zakharov;
    zak.MinVals = [-5, -5];
    zak.MaxVals = [10, 10];
    
    %Three Hump Camel
    thc = struct();
    thc.Name = 'ThreeHumpCamel';
    thc.Func = @ThreeHumpCamel;
    thc.MinVals = [-5, -5];
    thc.MaxVals = [5, 5];
    
    %Michalewicz Function
    mich = struct();
    mich.Name = 'Michalewicz';
    mich.Func = @Michalewicz;
    mich.MinVals = [0, 0];
    mich.MaxVals = [pi, pi];
    
    %Bukin Function No 6
    buk = struct();
    buk.Name = 'Bukin_6';
    buk.Func = @Bukin_6;
    buk.MinVals = [-15, -3];
    buk.MaxVals = [-5, 3];
    
    %Ackley Function 
    ack = struct();
    ack.Name = 'Ackley';
    ack.Func = @Ackley;
    ack.MinVals = [-32.77, -32.77];
    ack.MaxVals = [32.77, 32.77];
    
    %Levy Function No 13
    levy = struct();
    levy.Name = 'Levy_13';
    levy.Func = @Levi_13;
    levy.MinVals = [-10, -10];
    levy.MaxVals = [10, 10];
    
    %Eggholder Function
    egg = struct();
    egg.Name = 'Eggholder';
    egg.Func = @Eggholder;
    egg.MinVals = [-512, -512];
    egg.MaxVals = [512, 512];
    
    %Shekel Function
    shekel = struct();
    shekel.Name = 'Shekel';
    shekel.Func = @Shekel;
    shekel.MinVals = zeros(1,4);
    shekel.MaxVals = 10*ones(1,4);
    
    %Hartmann 6D Function
    hart6d = struct();
    hart6d.Name = 'Hartmann6D';
    hart6d.Func = @Hartmann6D;
    hart6d.MinVals = zeros(1,6);
    hart6d.MaxVals = ones(1,6);
    
    %Powell 12D Function
    powell10d = struct();
    powell10d.Name = 'Powell12D';
    powell10d.Func = @Powell;
    powell10d.MinVals = -4*ones(1, 12);
    powell10d.MaxVals = 5*ones(1, 12);
    
    %Ackley 5D Function
    ackley10d = struct();
    ackley10d.Name = 'Ackley5D';
    ackley10d.Func = @Ackley;
    ackley10d.MinVals = -32.77*ones(1, 5);
    ackley10d.MaxVals = 32.77*ones(1, 5);
    
    %Rosenbrock 20D Function
    rosen20d = struct();
    rosen20d.Name = 'Rosenbrock50D';
    rosen20d.Func = @RosenbrockNd;
    rosen20d.MinVals = repmat(-2.048, 1, 50);
    rosen20d.MaxVals = repmat(2.048, 1, 50);
    
    objFuns = {
        forr,
        sphere, 
        booth, 
        rosenbrock,
        beale,
        gp,
        zak,
        thc,
        mich,
        buk,
        ack,
        levy, 
        egg,
        shekel,
        hart6d,
        powell10d,
        ackley10d,
        rosen20d
        };
    
    


end