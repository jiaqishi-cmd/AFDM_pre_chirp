function metrics = evaluate_phi_metrics(H1, H2, delta)
%EVALUATE_PHI_METRICS Evaluate Phi(delta)=[H1*delta,H2*delta].

    phi1 = H1 * delta(:);
    phi2 = H2 * delta(:);
    Phi = [phi1, phi2];
    s = svd(Phi);

    metrics.rank = rank(Phi);
    metrics.sigma_min = min(s);
    metrics.sigma_max = max(s);
    metrics.sigma_ratio = metrics.sigma_min / max(metrics.sigma_max, eps);
    metrics.col_corr = abs(phi1' * phi2) / max(norm(phi1) * norm(phi2), eps);
    metrics.phi1 = phi1;
    metrics.phi2 = phi2;
end
