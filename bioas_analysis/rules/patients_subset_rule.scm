;; Patients Gene expression Subset introduction rule 
;; e.g 
;; EvaluationLink <TV>
;;   LazyExecutionOutputLink
;;     SchemaNode "make-overexpression-predicate-for-gene"  
;;     Gene
;;   patient
;; |-
;; SUbset <TV>
;;   Set patient
;;   SatisfyingSetScope
;;      Variable "$patient"
;;      EvaluationLink 
;;          LazyExecutionOutputLink
;;              SchemaNode "make-overexpression-predicate-for-gene"
;;              Gene
;;              Variable "$patient"

(define gene-expression-subset-rule
  (let* ((A (Variable "$A"))
         (B (Variable "$B"))
         (S (Variable "$S")))
    (Bind
      (VariableSet
        (TypedVariable A (Type "GeneNode"))
        (TypedVariable B (Type "PatientNode"))
        (TypedVariable S (Type "SchemaNode")))
      (Present
        (EvaluationLink
          (LazyExecutionOutputLink
            S
            A)
          B))
      (ExecutionOutput
        (GroundedSchema "scm: generate-subset")
        (List
          ;; conclusion
          (Subset 
            (Set B)
            (SatisfyingSetScope 
              (Variable "$patient")
              (EvaluationLink
                (LazyExecutionOutputLink
                  S
                  A)
                (Variable "$patient"))))
          ;; premises
          (EvaluationLink
            (LazyExecutionOutputLink
              S
              A)
            B))))))

(define patient-data-subset-rule
  (let* ((P (Variable "$P"))
         (Pr (Variable "$Pr"))
         (C (Variable "$C")))
    (Bind
      (VariableSet
        (TypedVariable P (Type "PatientNode"))
        (TypedVariable Pr (Type "PredicateNode"))
        (TypedVariable C (Type "ConceptNode")))
      (Present
        (EvaluationLink
          Pr
          (List
            P
            C)))
      (ExecutionOutput
        (GroundedSchema "scm: generate-subset")
        (List
          ;; conclusion
          (Subset 
            (Set P)
            (SatisfyingSetScope 
              (Variable "$patient")
              (EvaluationLink
                Pr
                (List
                  (Variable "$patient")
                  C))))
          ;; premises
          (EvaluationLink
            Pr
            (List
              P
              C)))))))

(define patient-data-boolean-subset-rule
  (let* ((P (Variable "$P"))
         (Pr (Variable "$Pr")))
    (Bind
      (VariableSet
        (TypedVariable P (Type "PatientNode"))
        (TypedVariable Pr (Type "PredicateNode")))
      (Present
        (EvaluationLink
          Pr
          P))
      (ExecutionOutput
        (GroundedSchema "scm: generate-subset")
        (List
          ;; conclusion
          (Subset 
            (Set P)
            (SatisfyingSetScope 
              (Variable "$patient")
              (EvaluationLink
                Pr
                (Variable "$patient"))))
          ;; premises
          (EvaluationLink
            Pr
            P))))))

(define (generate-subset conclusion . premises)
  (if (= (length premises) 1)
    (let* ((Ss conclusion)
          (eval (car premises))
          (st (cog-mean eval))
          (conf (cog-confidence eval)))
      (if (> conf 0) (cog-merge-hi-conf-tv! Ss (stv st conf))))))

; Name the rule
(define gene-expression-subset-rule-name
  (DefinedSchemaNode "gene-expression-subset-rule"))
(DefineLink gene-expression-subset-rule-name
  gene-expression-subset-rule)

(define patient-data-subset-rule-name
  (DefinedSchemaNode "patient-data-subset-rule"))
(DefineLink patient-data-subset-rule-name
  patient-data-subset-rule)

(define patient-data-boolean-subset-rule-name
  (DefinedSchemaNode "patient-data-boolean-subset-rule"))
(DefineLink patient-data-boolean-subset-rule-name
  patient-data-boolean-subset-rule)
