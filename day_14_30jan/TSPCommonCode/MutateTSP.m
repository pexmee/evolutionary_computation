function [Tour] = MutateTSP(Tour, MutationMethod)
% Mutate a closed TSP tour
NumCities = length(Tour) - 1;
if MutationMethod == 1
    % inversion
    StartNdx = randi([2, NumCities-1]);
    EndNdx = randi([StartNdx+1, NumCities]);
    Indices = [1:StartNdx-1, EndNdx:-1:StartNdx, EndNdx+1:NumCities+1];
    Tour = Tour(Indices);
elseif MutationMethod == 2
    % insertion
    FromPos = randi([2, NumCities-1]);
    while true
        ToPos = randi([2, NumCities-1]);
        if ToPos ~= FromPos, break, end
    end
    if ToPos < FromPos
        Indices = [1:ToPos-1, FromPos, ToPos:FromPos-1, FromPos+1:NumCities+1];
    else
        Indices = [1:FromPos-1, FromPos+1:ToPos-1, FromPos, ToPos:NumCities+1];
    end
    Tour = Tour(Indices);
elseif MutationMethod == 3
    % reciprocal exchange
    FromPos = randi([2, NumCities-1]);
    while true
        ToPos = randi([2, NumCities-1]);
        if ToPos ~= FromPos, break, end
    end
    if ToPos < FromPos
        Indices = [1:ToPos-1, FromPos, ToPos+1:FromPos-1, ToPos, FromPos+1:NumCities+1];
    else
        Indices = [1:FromPos-1, ToPos, FromPos+1:ToPos-1, FromPos, ToPos+1:NumCities+1];
    end
    Tour = Tour(Indices);
elseif MutationMethod == 4
    % opposition
    OppositeTour = zeros(NumCities, 1);
    for k = 1 : NumCities
        if mod(k, 2) == 1
            m = (k + 1) / 2;
        else
            m = NumCities + 1 - k / 2;
        end
        OppositeTour(k) = Tour(m);
    end
    Tour = [OppositeTour; OppositeTour(1)];
end
return