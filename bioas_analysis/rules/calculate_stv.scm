;; Calculate stv 

(define calculate_stv
  (let* ((Pt (Variable "$Pt"))
         (Ppty (Variable "$Ppty")))
    (Bind
      (VariableSet
        (TypedVariable Pt (Type "ConceptNode"))
        (TypedVariable Ppty (Type "SatisfyingSetScopeLink")))
      (Present
        (Subset
            (Set Pt)
            Ppty))
      (ExecutionOutput
              (GroundedSchema "scm: set-stv")
              (List
                (Subset (Set Pt) Ppty)
                (Set Pt)
                Ppty)))))

(define (set-stv Ss elem1 elem2)
  (let* ((Sstv (cog-tv Ss))
        (S1 (/ 1 total_patients))
        (conf (count->confidence total_patients))
        (res (cog-outgoing-set elem2))
        (S2 (/ (length (cog-outgoing-set (cog-execute! (Get res)))) total_patients)))
    (cog-merge-hi-conf-tv!
      (Subset 
        (cog-merge-hi-conf-tv! elem1 (stv S1 conf))
        (cog-merge-hi-conf-tv! elem2 (stv S2 conf))) 
      Sstv)))
