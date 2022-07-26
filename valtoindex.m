function index = valtoindex(val,numval,first,last)
index = round(((numval-1)/(last-first))*(val-first)+1);
end