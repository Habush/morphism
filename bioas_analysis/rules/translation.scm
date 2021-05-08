;; Crisp rules about translating a link into another link

;; Helpers
(define ConceptT (TypeInh "ConceptNode"))
(define GeneT (Type "GeneNode"))

(define-public (gen-present-link-translation-rule LINK-TYPE-1 LINK-TYPE-2 VAR-TYPE)
  (let* ((X (Variable "$X"))
         (Y (Variable "$Y"))
         (XY-1 (LINK-TYPE-1 X Y))
         (XY-2 (LINK-TYPE-2 X Y)))
    (Bind
      (VariableList 
        (TypedVariable X VAR-TYPE)
        (TypedVariable Y VAR-TYPE))
      (Present
        XY-1)
      XY-2)))

(define-public (inheritance->subset)
    (cog-logger-info "Inheritance->Subset")
    (gen-present-link-translation-rule InheritanceLink SubsetLink ConceptT)
    (cog-logger-info "Done!"))
  
