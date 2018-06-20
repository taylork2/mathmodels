medianIncome = VarName5;
states = TablewithrowheadingsincolumnAandcolumnheadingsinrows5to7;
natMedian = 5.1849*10^4;

size=length(medianIncome);
demGov = 1:19;
demSen = 1:17;
repGov = 1:31;
repSen = 1:20;

govCount = 1;
senCount = 1;
repGovCount = 1;
repSenCount = 1;

for i=1:size
    state = char(states(i));
    switch state
        case 'Alaska'
            continue
        case 'California'
        case 'Colorado'
        case 'Connecticut'
        case 'Delaware'
        case 'Hawaii'
        case 'Kentucky'
        case 'Minnesota'
        case 'Missouri'
        case 'Montana'
        case 'New Hampshire'
        case 'New York'
        case 'Oregon'
        case 'Pennsylvania'
        case 'Rhode Island'
        case 'Vermont'
        case 'Virginia'
        case 'Washington'
        case 'West Virginia'
        case 'District of Columbia'
        otherwise
            repGov(repGovCount) = medianIncome(i);
            repGovCount = repGovCount + 1;
            continue
    end
    demGov(govCount) = medianIncome(i);
    govCount = govCount + 1;
end

for i=1:size
    state = char(states(i));
    switch state
        case 'Colorado'
            continue
        case 'Florida'
            continue
        case 'Illinois'
            continue
        case 'Indiana'
            continue
        case 'Maine'
            continue
        case 'Missouri'
            continue
        case 'Montana'
            continue
        case 'Nevada'
            continue
        case 'New Hampshire'
            continue 
        case 'North Dakota'
            continue
        case 'Ohio'
            continue
        case 'Pennsylvania'
            continue
        case 'West Virginia'
            continue
        case 'Wisconsin'
            continue
        case 'California'
        case 'Connecticut'
        case 'Delaware'
        case 'Hawaii'
        case 'Maryland'
        case 'Minnesota'
        case 'Massachusetts'
        case 'Michigan'
        case 'New Jersey'
        case 'New Mexico'
        case 'New York'
        case 'Oregon'
        case 'Rhode Island'
        case 'Vermont'
        case 'Virginia'
        case 'Washington'
        case 'District of Columbia'
        otherwise
            repSen(repSenCount) = medianIncome(i);
            repSenCount = repSenCount + 1;
            continue
    end
    demSen(senCount) = medianIncome(i);
    senCount = senCount + 1;
end

govSize = length(demGov);
senSize = length(demSen);
repGovSize = length(repGov);
repSenSize = length(repSen);
govMu = mean(demGov);
senMu = mean(demSen);
repSenMu = mean(repSen);
repGovMu = mean(repGov);

repSenSigma=(sum((repSen-repSenMu.*ones(1, repSenSize)).^2)/(repSenSize))^(1/2);
SE1 = repSenSigma*sqrt(repSenSize);
z1 = (repSenMu - natMedian)/SE1

senSigma=(sum((demSen-senMu.*ones(1, senSize)).^2)/(senSize))^(1/2);
SE2 = senSigma*sqrt(senSize);
z2 = (senMu - natMedian)/SE2

repGovSigma=(sum((repGov-repGovMu.*ones(1, repGovSize)).^2)/(repGovSize))^(1/2);
SE3 = repGovSigma*sqrt(repGovSize);
z3 = (repGovMu - natMedian)/SE3

govSigma=(sum((demGov-govMu.*ones(1, govSize)).^2)/(govSize))^(1/2);
SE4 = govSigma*sqrt(govSize);
z4 = (govMu - natMedian)/SE4

    
