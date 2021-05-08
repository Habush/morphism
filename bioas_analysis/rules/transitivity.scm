
;; Crisp rules about transitivity of predicates or inheritance links

;; Helpers
(define ConceptT (TypeInh "ConceptNode"))
(define GeneT (Type "GeneNode"))

(define-public (gen-present-link-transitivity-rule LINK-TYPE VAR-TYPE)
  (let* ((X (Variable "$X"))
         (Y (Variable "$Y"))
         (Z (Variable "$Z"))
         (YZ (LINK-TYPE Y Z))
         (XZ (LINK-TYPE X Z)))
    (Bind
      (VariableList
        (TypedVariable X VAR-TYPE)
        (TypedVariable Y VAR-TYPE)
        (TypedVariable Z VAR-TYPE))
      (And
        (Present
          XY
          YZ)
        (Not (Identical X Z)))
      XZ)))

(define-public (gen-present-mixed-link-transitivity-rule LINK-TYPE-1 LINK-TYPE-2
                                                  X-TYPE Y-TYPE Z-TYPE)
  (let* ((X (Variable "$X"))
         (Y (Variable "$Y"))
         (Z (Variable "$Z"))
         (XY (LINK-TYPE-1 X Y))
         (YZ (LINK-TYPE-2 Y Z))
         (XZ (LINK-TYPE-1 X Z)))
    (Bind
      (VariableList
        (TypedVariable X X-TYPE)
        (TypedVariable Y Y-TYPE)
        (TypedVariable Z Z-TYPE))
      (And
        (Present
          XY
          YZ)
        (Not (Identical X Z)))
      XZ)))

(define-public (gen-present-link-transitivity)
    (cog-logger-info "Running gen-present-link-transitivity")
    (gen-present-link-transitivity-rule SubsetLink ConceptT)
    (cog-logger-info "Done!")))

(define-public (gen-present-mixed-link-transitivity)
    (gen-present-mixed-link-transitivity-rule MemberLink
                                            SubsetLink
                                            GeneT ConceptT ConceptT)
    (cog-logger-info "Done!"))