namespace Calc_Energy {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;

    operation QPE(t : Int, U : (Qubit[], Int) => Unit is Adj + Ctl, q : Qubit[], idx : Int) : Double {
        use q_mas = Qubit[t];
 for i in 0..t - 1 {
            H(q_mas[i]);
        }
        
        for i in 0..t - 1 {
            let power = 2^i;
            for j in 1..power {
                Controlled U([q_mas[i]], (q, idx));
            }
        }

        SwapReverseRegister(q_mas);    
        Adjoint ApplyQFT(q_mas);

        let result = MeasureInteger(q_mas);
        ResetAll(q_mas);
        return IntAsDouble(result) / (2.0^IntAsDouble(t));
        
    }

    function H2Coeffs(idx : Int) : Double[] {
        return [
            [0.2,  -0.4803,  0.8919, -0.8919,  0.1224,  0.5328, 0.5328],
            [0.4,  -0.8952,  0.6984, -0.6984,  0.0748,  0.4330, 0.4330],
            [0.6,  -1.0513,  0.5214, -0.5214,  0.0339,  0.3109, 0.3109],
            [0.735,-1.0524,  0.3979, -0.3979, -0.0112,  0.1809, 0.1809],
            [0.8,  -1.0400,  0.3603, -0.3603, -0.0329,  0.1509, 0.1509],
            [1.0,  -0.9584,  0.2622, -0.2622, -0.0785,  0.0979, 0.0979],
            [1.2,  -0.8766,  0.1945, -0.1945, -0.1053,  0.0679, 0.0679],
            [1.4,  -0.8062,  0.1551, -0.1551, -0.1229,  0.0469, 0.0469],
            [1.6,  -0.7468,  0.1235, -0.1235, -0.1324,  0.0341, 0.0341],
            [1.8,  -0.6966,  0.0998, -0.0998, -0.1372,  0.0254, 0.0254],
            [2.0,  -0.6539,  0.0816, -0.0816, -0.1391,  0.0192, 0.0192],
            [2.5,  -0.5654,  0.0519, -0.0519, -0.1393,  0.0105, 0.0105],
            [3.0,  -0.5044,  0.0358, -0.0358, -0.1360,  0.0063, 0.0063],
            [4.0,  -0.4312,  0.0198, -0.0198, -0.1287,  0.0025, 0.0025]
        ][idx];
    }

    operation H2Evolution(q : Qubit[], idx : Int) : Unit is Adj + Ctl {
        let pauli = [
            [PauliI, PauliI],
            [PauliZ, PauliI], 
            [PauliI, PauliZ], 
            [PauliZ, PauliZ], 
            [PauliY, PauliY], 
            [PauliX, PauliX]
        ];

        for i in 0..5 {
            let coeff = H2Coeffs(idx)[i + 1];
            Exp(pauli[i], -2.0 * coeff, q);
        }
    }
    

    operation Calc(idx : Int) : Double {
        let t = 10;
        use q = Qubit[2];
        X(q[0]);
        X(q[1]);
        let phase = QPE(t, H2Evolution, q, idx);
        ResetAll(q);
        return -2.0 * PI() * phase / 2.0;
    }
}
       


