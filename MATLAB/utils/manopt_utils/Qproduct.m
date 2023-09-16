function QX = Qproduct(X, problem_data)
    if ~problem_data.use_marginalized
        QX = problem_data.Q * X;
    else
        QmainX = problem_data.Qmain * X;

        QtransX = problem_data.LeftOperatorRed * ...
         (problem_data.LtransCholRed' \ ...
            (problem_data.LtransCholRed \ ...
                (problem_data.LeftOperatorRed' * X)));

        QX = QmainX - QtransX;
    end

end

