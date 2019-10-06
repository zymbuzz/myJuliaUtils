function dataTransf(x::AbstractVecOrMat, transfId::Vector{Int})
    y = similar(x)
    n_loss = 0
    for i = 1:size(x, 2)
        (y[:,i], n_loss_temp) = dataTransf(x[:,i], transfId[i])
        n_loss = max(n_loss, n_loss_temp)
    end
    return (y = y, n_loss = n_loss)
end

function dataTransf(x::AbstractVecOrMat, transfId::Int)
#   Categories:
#   ------------------------------------------------------------------------
#   1. Levels (Y_{t})
#   2. First Difference (Y_{t}-Y_{t-1})
#   3. Seasonal difference in quarterly data ((Y_{t})-(Y_{t-4}))
#   4. Seasonal difference in monthly data ((Y_{t})-(Y_{t-12}))
#   5. Log-levels (Ln(Y_{t}))
#   6. Log-First Difference (Ln(Y_{t})-Ln(Y_{t-1}))
#   7. Seasonal log difference in quarterly data (Ln(Y_{t})-Ln(Y_{t-4}))
#   8. Seasonal log difference in monthly data (Ln(Y_{t})-Ln(Y_{t-12}))
#   9. Detrending Ln(Y_{t}) by HP filter using quarterly data
#   10. Detrending Ln(Y_{t}) by HP filter using monthly data
#   11. Seasonal difference of a detrended Ln(Y_{t}) by HP filter using quarterly data
#   12. Seasonal difference of a detrended Ln(Y_{t}) by HP filter using monthly data
#   13. Seasonal difference of a detrended Ln(Y_{t}) by removing a linear trend
#   14. Detrended Y_{t} by removing a linear trend (for already SA data)
#   15. Detrended Ln(Y_{t}) by removing a linear trend (for already SA data)
#  =========================================================================
    y = similar(x)

    if transfId == 1
        y = copy(x)
        n_loss = 0
    elseif transfId == 2
        y[2:end,:] = x[2:end,:] - x[1:end - 1,:]
        n_loss = 1
    elseif transfId == 3
        y[5:end,:] = x[5:end,:] - x[1:end - 4,:]
        n_loss = 4
    elseif transfId == 4
        y[13:end,:] = x[13:end,:] - x[1:end - 12,:]
        n_loss = 12
    elseif transfId == 5
        y = log.(x)
        n_loss = 0
    elseif transfId == 6
        y[2:end,:] = log.(x[2:end,:]) - log.(x[1:end - 1,:])
        n_loss = 1
    elseif transfId == 7
        y[5:end,:] = log.(x[5:end,:]) - log.(x[1:end - 4,:])
        n_loss = 4
    elseif transfId == 8
        y[13:end,:] = log.(x[13:end,:]) - log.(x[1:end - 12,:])
        n_loss = 12
        # case 7
        #     yt=hpfilter(log(x),1600);
        #     y=log(x)-yt;
        #     n_loss=0;
        # case 8
        #     yt=hpfilter(log(x),14400);
        #     y=log(x)-yt;
        #     n_loss=0;
        # case 9
        #     yt=hpfilter(log(x),1600);
        #     yt=log(x)-yt;
        #     y(5:end,:)=yt(5:end,:)-yt(1:end-4,:);
        #     n_loss=4;
        # case 10
        #     yt=hpfilter(log(x),14400);
        #     yt=log(x)-yt;
        #     y(13:end,:)=yt(13:end,:)-yt(1:end-12,:);            
        #     n_loss=12;
    elseif transfId == 13
        yt = log.(x)
        yt = yt[13:end,:] - yt[1:end - 12,:]           
        y[14:end,:] = yt[2:end,:] - yt[1:end - 1,:]
        n_loss = 13
    elseif transfId == 14
        y = detrend(x, power=1)
        n_loss = 0 
    elseif transfId == 15
        y = detrend(log.(x), power=1)
        n_loss = 0 
    else
        error("Please specify a valid transformation index")
    end

    return (y = y, n_loss = n_loss)
end

function detrend(y::AbstractVecOrMat;power::Int = 1)
    T = size(y, 1)
    X = similar(y, T, power + 1)
    X[:,1] = ones(T)
    if power > 0
        for i = 1:power
            X[:,i + 1] = collect(1:T).^i
        end
    end
    b = X \ y
    return y - X * b
end