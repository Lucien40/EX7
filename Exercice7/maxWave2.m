function [A,u,tmax] = maxWave2(x,t,f,peak,fitrange,lim)
index = (1:length(t))';
tmax = zeros(length(x),1);
A = zeros(length(x),1);
for i = 1:length(x)
    fi = f(:,i);
    lb = 1;
    for j = 1:peak
        if rem(j,2) ~= 0
            ilim = index(index >= lb & fi > lim);
        else
            ilim = index(index >= lb & fi < -lim);
        end
        lb = ilim(1);
    end
    if rem(j,2) ~= 0
        ilim = index(index >= lb & fi <= lim);
    else
        ilim = index(index >= lb & fi >= -lim);
    end
    rb = ilim(1) - 1;
    [~,imax] = max((-1)^(peak + 1) * fi(lb:rb));
    imean = floor(mean(imax));
    lfitb = max(lb,lb - 1 + imean - fitrange);
    rfitb = min(rb,lb - 1 + imean + fitrange);
    p = polyfit(t(lfitb:rfitb),fi(lfitb:rfitb),2);
    tmax(i) = -p(2)/(2*p(1));
    A(i) = (-1)^(peak + 1) * polyval(p,tmax(i));
end
dx = diff(x);
dtmax = diff(tmax);
u = dx./dtmax;
end

