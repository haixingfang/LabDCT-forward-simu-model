function multi = find_multiplicity(pos,sg)

for i=1:sg.nouniq
    lp(i,:) = pos*sg.symop(i).rot+sg.symop(i).trans;
end

lpu(1,:) = lp(1,:);

for i=2:sg.nouniq
    nuniq = size(lpu,1);
    for j=1:nuniq
        t1 = lp(i,1)-lpu(j,1);
        t2 = lp(i,2)-lpu(j,2);
        t3 = lp(i,3)-lpu(j,3);
        if mod(t1,1)< 0.0001 & mod(t2,1)< 0.0001 & mod(t3,1)< 0.0001
            break
        else
            if j == nuniq
                lpu(nuniq+1,:) = lp(i,:);
            end
        end
    end
end
uniq = size(lpu,1);
multi = sg.nouniq/uniq;
