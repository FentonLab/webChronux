%get sliding variance
%winLen = half window (sample -/+ winLen)
sv = rvar(data,winLen,length(data));