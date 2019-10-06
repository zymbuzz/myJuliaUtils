
function runtvpvarMLTVflat(data, REPS::Int64, KEEP::Int64, KEEPgap::Int64, it_print::Int64, CONST::Int64, P::Int64, scaleb0var::Real, scaleQ0::Real)
    # data
    # REPS::Int64
    # KEEP::Int64
    # KEEPgap::Int64
    # it_print::Int64
    # CONST::Int64
    # P::Int64
    # scaleb0var::Real
    # scaleQ0::Real

    BURN = REPS - KEEP * KEEPgap

    y, x = preparexy(data, P, CONST)
    T, n = size(y)
    NP = size(x, 2)
    K = n * NP # number of elements in the state vector

    # uninformative priors          
    b0 = zeros(K, 1)
    b0V = scaleb0var .* eye(K)

    # Q is the covariance of B(t)
    # Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
    Q0 = scaleQ0 .* eye(K)
    # Q0 = zeros(K, K);
    Q0V = K + 1
    QVdf = T + Q0V

    # Sigma is the covariance of the VAR covariance, SIGMA
    # Sigma ~ IW(I,n+1)
    # sigma0 = zeros(n, n);
    sigma0V = n + 1
    sigmaVdf = T + sigma0V

    y = y' 
    y = convert(Array{Float64,2}, y);
    Z = kron(x, eye(n));

    # initialising:
    Qdraw = 0.01 .* eye(K);
    sigmadraw = 0.1 .* eye(n);
    bdraw = Array{Float64}(undef, K, T)

    # Storage matrices for posteriors and stuff
    bmean = zeros(K, T);
    Qmean = zeros(K, K);
    sigmamean = zeros(n, n);

    bmean = zeros(BigFloat, K, T);
    Qmean = zeros(BigFloat, K, K);
    sigmamean = zeros(BigFloat, n, n);

    iter_store = 1;

    no2bestored = min(50, REPS);
    bstore = zeros(K * T, no2bestored);
    Qstore = zeros(K * K, no2bestored);
    sigmastore = zeros(n * n, no2bestored);
    btstore = zeros(K * T, no2bestored);
    Vtstore = zeros(K^2 * T, no2bestored);

    for iter = 1:REPS

        # step 1, drawing b
        bdraw, bt, Vt = carterkohnMLTV(y, Z, sigmadraw, Qdraw, K, n, T, b0, b0V);

        # step 2, drawing Q
        Btemp = diff(bdraw, dims=2);
        sse_2Q = Btemp * Btemp';

        # Qdraw=rand(InverseWishart( T + Q0V,sse_2Q + Q0))

        Qinv = invss(sse_2Q + Q0)
        Qinvdraw = wish(Qinv, QVdf)
        Qdraw = inv(Qinvdraw)

        # step 3, drawing sigma
        sse_2S = zeros(n, n);
        yhat = Array{Float64}(undef, n)
        for i = 1:T
            yhat = y[:, i] - Z[(i - 1) * n + 1:i * n,:] * bdraw[:, i];
            sse_2S += yhat * yhat';
        end   

        # sigmadraw=rand(InverseWishart( T + sigma0V,sse_2S + sigma0));

        Sigmainv = invss(sse_2S);
        Sigmainvdraw = wish(Sigmainv, sigmaVdf);
        sigmadraw = inv(Sigmainvdraw);
        
        if mod(iter, it_print) == 0
            println("Iter: $iter / $REPS ")
            # println(bdraw)
            # println(Qdraw)
            # println(sigmadraw)

            # println(cov(Btemp'))

            # println(bmean)
            # println(Qmean)
            # println(sigmamean)
        end

        # Saving draws after burn
        if iter > BURN && iter_store == KEEPgap;
            iter_store = 1;
            bmean += bdraw;
            Qmean += Qdraw;
            sigmamean += sigmadraw;
        elseif iter > BURN
            iter_store += 1;
        end

        # show(vec(bdraw))
        # sleep(20)
        # show(vec(Qdraw))
        # sleep(20)
        # show(vec(sigmadraw))
        # sleep(20)

        if iter <= no2bestored
            bstore[:,iter] = vec(bdraw);
            Qstore[:,iter] = vec(Qdraw);
            sigmastore[:,iter] = vec(sigmadraw);
            btstore[:,iter] = vec(bt);
            Vtstore[:,iter] = vec(Vt);

        end

    end
    # computing averages
    bmean = bmean ./ KEEP;
    Qmean = Qmean ./ KEEP;
    sigmamean = sigmamean ./ KEEP; 

    # writedlm('StoredDraws/bstore2',bstore,',')
    # writedlm('StoredDraws/Qstore2',Qstore,',')
    # writedlm('StoredDraws/sigmastore2',sigmastore,',')

    # return bmean, Qmean, sigmamean
    return bmean, Qmean, sigmamean, bstore, Qstore, sigmastore, btstore, Vtstore
end
