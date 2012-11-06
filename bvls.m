function [x, iter] = bvls(A, b, l, u)

iterMax = 500;
epsilon = 1e-12;

SETF = 0;
SETL = 1;
SETU = 2;

m = size(A,1);
n = size(A,2);
check1Failed = false;
iAlphaMin = -1;
normb = norm(b);

% assume cold start
for i=1:n
    set(i) = SETL;
end
x = l;

iter = 0;

while (iter < iterMax)
    if (norm(A*x-b) <= normb*epsilon)
        break;
    end
    
    if (~check1Failed)
        w = A'*(b-A*x);
    end
    
    % check 2
    if (iAlphaMin >= 0)
        w(iAlphaMin) = 0;
    end
    iAlphaMin = -1;
    
    isAllInF = true;
    isAllwLNeg = true;
    isAllwUPos = true;
    for i=1:n
        if ((set(i) == SETL) && (isAllInF || isAllwLNeg))
            isAllInF = false;
            if (w(i) > 0)
                isAllwLNeg = false;
            end
        elseif ((set(i) == SETU) && (isAllInF || isAllwUPos))
            isAllInF = false;
            if (w(i) < 0)
                isAllwUPos = false;
            end
        end
    end
    if (isAllInF || (isAllwLNeg && isAllwUPos))
        break;
    end
    
    tstar = -1;
    for i=1:n
        if (set(i) ~= SETF)            
            if (set(i) == SETL)
                sw = w(i);
            elseif (set(i) == SETU)
                sw = -w(i);
            end
            if ((tstar < 0) || (sw > swMax))
                tstar = i;
                setTstar = set(tstar);
                swMax = sw;
            end
        end
    end
    
    set(tstar) = SETF;
    
    while (iter < iterMax)
        Ap = [];
        bp = b;
        for i=1:n
            if (set(i) ~= SETF)
                bp = bp - A(:,i)*x(i);
            else
                Ap = [Ap A(:,i)];
            end
        end
        if isempty(Ap)
            break;
        end
        z = pinv(Ap)*bp;
        
        % check 1
        check1Failed = false;
%         if (rank(Ap) < size(Ap,2))
%             set(tstar) = setTstar;
%             w(tstar) = 0;
%             tstar = -1;
%             check1Failed = true;
%         else
        i = 1;
        for j=1:tstar
            if (set(j) == SETF)
                if (j == tstar)
                    if ( ((setTstar == SETL) && (z(i) < l(j))) ||...
                         ((setTstar == SETU) && (z(i) > u(j))) )
                        set(tstar) = setTstar;
                        w(tstar) = 0;
                        tstar = -1;
                        check1Failed = true;
                    end
                end
                i = i+1;
            end
        end
%         end
        if (check1Failed)
            break;
        end
        
        cond = true;
        i = 1;
        for j=1:n
            if (set(j) == SETF)
                if ((z(i) <= l(j)) || (z(i) >= u(j)))
                    cond = false;
                    break;
                end
                i = i+1;
            end
        end
        if cond
            i = 1;
            for j=1:n
                if (set(j) == SETF)
                    x(j) = z(i);
                    i = i+1;
                end
            end
            iter = iter+1;
            break;
        end        
        alphaMin = -1;
        iAlphaMin = -1;
        sAlphaMin = 0;
        i = 1;
        for j=1:n
            if (set(j) == SETF)
                if ((z(i) <= l(j)) || (z(i) >= u(j)))
                    if z(i) >= u(j)
                        alpha = (u(j) - x(j))/(z(i) - x(j));
                    else
                        alpha = (l(j) - x(j))/(z(i) - x(j));
                    end
                    if ((alphaMin < 0) || (alpha < alphaMin))
                        alphaMin = alpha;
                        iAlphaMin = j;
                        if z(i) < l(j)
                            sAlphaMin = -1;
                        else
                            sAlphaMin = 1;
                        end
                    end
                end
                i = i+1;
            end
        end
        i = 1;
        for j=1:n
            if (set(j) == SETF)
                x(j) = x(j) + alphaMin*(z(i)-x(j));
                i = i+1;
            end
        end
        if sAlphaMin > 0
            x(iAlphaMin) = u(iAlphaMin);
        else
            x(iAlphaMin) = l(iAlphaMin);
        end
        for i=1:n
            if (x(i) <= l(i))
                set(i) = SETL;
            elseif (x(i) >= u(i))
                set(i) = SETU;
            end
        end
        iter = iter+1;
    end
end
