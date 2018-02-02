;;; qmd-progress --- Fortran beautification settings.
;;;
;;; This file contains formatting settings for the Emacs f90-mode. Add
;;; it to your .emacs file.

(add-hook 'f90-mode-hook
          ;; Extra indentation applied to DO blocks (default 3).
          (setq-local f90-do-indent 3)

          ;; Extra indentation within if/select/where/forall blocks (default 3).
          (setq-local f90-if-indent 3)

          ;; Extra indentation within type/enum/interface/block-data
          ;; blocks (default 3).
          (setq-local f90-type-indent 3)

          ;; Extra indentation within
          ;; program/module/subroutine/function blocks (default 2).
          (setq-local f90-program-indent 2)

          ;; Extra indentation within associate blocks (default 2).
          (setq-local f90-associate-indent 2)

          ;; Extra indentation within critical/block blocks (default
          ;; 2).
          (setq-local f90-critical-indent 2)

          ;; Extra indentation applied to continuation lines (default
          ;; 5).
          (setq-local f90-continuation-indent 5)

          ;; String inserted by function M-x f90-comment-region at
          ;; start of each line in region (default "!!!$").
          (setq-local f90-comment-region "!!!$")

          ;; Regexp determining the type of comment to be intended
          ;; like code (default "!").
          (setq-local f90-indented-comment-re "!")

          ;; Regexp of comment-like directive like "!HPF\\$", not to
          ;; be indented (default "!hpf\\$").
          (setq-local f90-directive-comment-re "!hpf\\$")

          ;; Regexp holding list of delimiters at which lines may be
          ;; broken (default "[-+*/><=,% \t]").
          (setq-local f90-break-delimiters "[-+*/><=,% \t]")

          ;; Non-nil causes ‘f90-do-auto-fill’ to break lines before
          ;; delimiters (default t).
          (setq-local f90-break-before-delimiters t)

          ;; Automatic insertion of ‘&’ at beginning of continuation
          ;; lines (default t).
          (setq-local f90-beginning-ampersand t)

          ;; From an END statement, check and fill the end using
          ;; matching block start. Allowed values are ‘blink’,
          ;; ‘no-blink’, and nil, which determine whether to blink the
          ;; matching beginning (default ‘blink’).
          (setq-local f90-smart-end blink)

          ;; Automatic change of case of keywords (default nil). The
          ;; possibilities are ‘downcase-word’, ‘upcase-word’,
          ;; ‘capitalize-word’.
          (setq-local f90-auto-keyword-case nil)

          ;; Do not left-justify line numbers (default nil).
          (setq-local f90-leave-line-no nil)
          )
