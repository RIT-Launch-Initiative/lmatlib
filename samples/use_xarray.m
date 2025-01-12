clear;

xarr = xarray(randi(100, 10, 50, 3), ...
    time = seconds(1:10), x = linspace(0, 10, 50), ...
    name = ["first", "second", "third"])
another_xarr = xarray(randi(100, 100, 50, 3), ...
    time = seconds(11:110), x = linspace(0, 10, 50), ...
    name = ["first", "second", "third"]) 
xarr = cat("time", xarr, another_xarr)
xarr = permute(xarr, ["name", "x", "time"])
xarr = xarr.range(x = [2 5]) % select by membership in range
xarr = xarr.pickt(time = [seconds(4), seconds(1)], x = [4 0.5]) % select by tolerance around value
xarr = xarr.pick(name = "second") % select by 
another_xarr = another_xarr.range(x = [2 5]).pickt(time = [seconds(20) seconds(2)]) % chain selections

