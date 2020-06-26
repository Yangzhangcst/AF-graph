function result = histRate(x) 

    if isnumeric(x) 
        x = x(:); 
        x = x(~isnan(x)); 
        xid = []; 
        else 
        [x,xid] = grp2idx(x); 
        x = x(~isnan(x)); 
    end 

    x = sort(x(:)); % ?? 
    m = length(x); 
    x1 = diff(x); % ??? 
    x1(end + 1) = 1; 
    x1 = find(x1); 
    CumFreq = x1/m; 
    value = x(x1); 
    x1 = [0; x1]; 
    Freq1 = diff(x1); 
    Freq2 = Freq1/m; 
    if nargout == 0 
    if isempty(xid) 
    fmt1 = '%11s %8s %6s %6s\n'; 
    fmt2 = ' %10d %8d %6.2f%% %6.2f%%\n'; 
    fprintf(1, fmt1, '??', '??', '??', '????'); 
    fprintf(1, fmt2, [value'; Freq1'; 100*Freq2'; 100*CumFreq']); 
    else 
    head = {'??', '??', '??(%)', '????(%)'}; 
    [head;xid,num2cell([Freq1, 100*Freq2, 100*CumFreq])] 
    end 
    else 
    if isempty(xid) result = [value Freq1 Freq2 CumFreq]; 
    else 
    result = [xid,num2cell([Freq1, Freq2, CumFreq])]; 
    end 
end 