# Documentation

See the [manual on readthedocs](https://meep.readthedocs.io/en/latest) for the latest documentation.

How to use magnetized plasma
(define-param w_p_meep 1); plasma angular frequency in meep units (do not forget to use 2*pi factor)
(define-param w_b_meep 1); cyclotron angular frequency in meep units (do not forget to use 2*pi factor)
(define-param v_c_meep 1); collision frequency in meep units
(define-param semicond_off1
   (make medium (epsilon 1.0 )
 	  (E-susceptibilities
		(make gyroelectric-susceptibility (omega_plas w_p_meep) (nu_cycl w_b_meep) (nu_col v_c_meep)   (sigma-offdiag w_p_meep 0 0) (b_dir 0 0 1)) 
	   )
	)
)
