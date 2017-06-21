%compute modulation index as shown in Tort 2010
%Dino Dvorak 2012 dino@indus3.net

function MI = GetMI(amp, phase, edges)

%get histogram for phases
h = zeros(1,length(edges)-1);
for hi = 1:length(edges) - 1
    k = phase >= edges(hi) & phase < edges(hi+1);
    if sum(k) > 0 %only if there are some values, otherwise keep zero
        h(hi) = mean(amp(k)); %mean of amplitudes at that phase bin
    end
end

%fix last value
k = find(phase == edges(end));
if ~isempty(k)
    h(end) = h(end) + mean(amp(k));
end

%convert to probability
h = h / sum(h);

%replace zeros by eps
k = h == 0;
h(k) = eps;

%calculate modulation index
Hp = -1 * sum(h .* log(h)); %entropy of the histogram
D = log(length(h)) - Hp;
MI = D / log(length(h));


