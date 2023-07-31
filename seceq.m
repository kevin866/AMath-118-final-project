function [root, num_iter] = seceq(n, i, b, d, tol)
    num_iter = 0;
    x = (b(1)*d(2)+d(1)*b(2))/(b(1)+b(2));
    if i > 1
        offset = sum(b(1:i-1)./(d(1:i-1)-d(2:i)));
        c = -((d(i)+d(i+1))*offset+(b(i)+b(i+1)));
        e = b(i)*d(i+1)+d(i)*b(i+1)+d(i)*d(i+1)*offset;
        x = (-c-sqrt(c^2-4*offset*e))/(2*offset);
    end
    b1 = b(1:i);
    b2 = b(i+1:n);
    d1 = d(1:i);
    d2 = d(i+1:n);
    while num_iter < 10
        phi = 1+sum(b1./(d1-x));
        phi_dif = sum(b1./((d1-x).^2));
        beta = -phi_dif*x^2;
        alpha = phi-beta/x;
        psi = sum(b2./(d2-x));
        psi_dif = sum(b2./((d2-x).^2));
        delta = x + psi/psi_dif;
        gamma = psi*(delta-x);
        q = beta-alpha*delta-gamma;
        root = (-q-sqrt(q^2+4*alpha*beta*delta))/(2*alpha);
        if abs((root-x)/root) <= tol
            break
        end
        x = root;
        num_iter = num_iter+1;
    end
end