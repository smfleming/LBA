function [v A b sv t0] = LBA_parse(model, pArray, Ncond)
% Parse parameter vector for fitting LBA
%
% SF 2012


%% Parse pArray vector and setup some transformations
j=1;
err = 0;
if model.v == 1
    v = repmat(pArray(j), Ncond, 1);
    j=j+1;
elseif model.v == Ncond
    for c = 1:model.v
        v(c,:) = pArray(j);
        j=j+1;
    end
else err = 1;
end

if model.A == 1
    A = repmat(pArray(j), Ncond, 1);
    j=j+1;
elseif model.A == Ncond
    for c = 1:model.A
        A(c,:) = pArray(j);
        j=j+1;
    end
else err = 1;
end

if model.b == 1
    b = repmat(pArray(j), Ncond, 1);
    j=j+1;
elseif model.b == Ncond
    for c = 1:model.b
        b(c,:) = pArray(j);
        j=j+1;
    end
else err = 1;
end

if model.sv == 1
    sv = repmat(pArray(j), Ncond, 1);
    j=j+1;
elseif model.sv == Ncond
    for c = 1:model.sv
        sv(c,:) = pArray(j);
        j=j+1;
    end
else err = 1;
end

if length(pArray) >= j
    if model.t0 == 1
        t0 = repmat(pArray(j), Ncond, 1);
        j=j+1;
    elseif model.t0 == Ncond
        for c = 1:model.t0
            t0(c,:) = pArray(j);
            j=j+1;
        end
    else err = 1;
    end
else
    t0 = repmat(200, Ncond, 1); % fix if not free param
    j=j+1;
end

if err == 1
    fprintf('\n\n\nBad input! model fields can only take on values of 1 or Ncond. See help LBA_mle.\n\n');
    return;
elseif length(pArray) > j-1
    fprintf('\n\n\nBad input! pArray does not match model specification. See help LBA_mle.\n\n');
    return;
end
