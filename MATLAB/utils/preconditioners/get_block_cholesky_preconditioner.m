function precon_function = get_block_cholesky_preconditioner(prob)
    rot_idxs = prob.all_R_idxs;
    Qrot = prob.Q(rot_idxs, rot_idxs);
    Lrot = chol(Qrot, 'lower');
    LrotT = Lrot';

    dist_idxs = prob.all_d_idxs;
    Qdist = prob.Q(dist_idxs, dist_idxs);
    QdistInv = sparse(diag(1./diag(Qdist)));
    Ldist = chol(Qdist, 'lower');

    Lmain = blkdiag(Lrot, Ldist);
    LmainT = Lmain';

    if prob.use_marginalized
        precon_function = @(x) LmainT \ (Lmain \ x);
    else

        % ignore last translation index, we will effectively be pinning it
        tran_idxs = [prob.all_t_idxs, prob.all_l_idxs];
        % tran_idxs = tran_idxs(2:end);
        tran_idxs = tran_idxs(1:end-1);

        Qtran = prob.Q(tran_idxs, tran_idxs);
        [Ltran, flag, Ptran] = chol(Qtran, 'lower');
        LtranT = Ltran';
        nmain = size(Lmain, 1);
        % precon_function = @(x) [LT \ (L \ x(1:Q_width-1, :)); zeros(1, size(x, 2))];
        precon_function = @(x) [LmainT \ (Lmain \ x(1:nmain, :));...
                                % zeros(1, size(x, 2));...
                                Ptran * (LtranT \ (Ltran \ (Ptran' * x(tran_idxs, :))));...
                                zeros(1, size(x, 2))...
                                ];
        % A = (P * (L2' \ (L2 \( P' * Qtran))));

        % precon_function = @(x) [LrotT \ (Lrot \ x(rot_idxs, :));...
        %                         QdistInv * x(dist_idxs, :);...
        %                         % zeros(1, size(x, 2));...
        %                         Ptran * (LtranT \ (Ltran \ (Ptran' * x(tran_idxs, :))));...
        %                         zeros(1, size(x, 2))...
        %                     ];

    end

end