;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 20:19:57 2008 (CEST)
;; Last-Updated: 2011-08-11
;;     Update #: 7

(in-package :cltmhmm)

(defun identify-region (line)
  (labels ((subchain (type line pos)
             (setq line (remove-if #'(lambda(x) (char= x #\tab)) line))
             ;;always suppose the second number equal or bigger than the first; 2 2 -> size = 1
             (let ((posend 0) (start 0) (end 0))
               (setq pos (position-if #'(lambda(x) (digit-char-p x)) line
                                          :start (+ pos (case type
                                                          (#\i 6)
                                                          ((#\o #\M) 7)
                                                          (t 0)))))
               (setq posend (position #\space line :start pos))
               (setq start (string-to-number (subseq line pos posend)))
               (setq pos (position-if #'(lambda(x) (digit-char-p x)) line :start posend))
               (setq posend (position-if-not #'(lambda(x) (digit-char-p x)) line :start pos))
               (setq end (string-to-number (subseq line pos posend)))
               (make-string (1+ (- end start)) :initial-element type))))
    (let ((pos 0))
      (setq line (substitute #\space #\tab line))
      (cond
        ((setq pos (search "inside" line)) (subchain #\i line pos))
        ((setq pos (search "TMhelix" line)) (subchain #\M line pos))
        ((setq pos (search "outside" line)) (subchain #\o line pos))
        (t "")))))

(defun translate-tmhmm-results (file &optional (command-symbol #\#))
  (with-open-file (stream file)
    (do ((predictions)
         (aux)
         (line (read-line stream nil nil) (read-line stream nil nil)))
        ((null line) (when (string/= "" aux) (push aux predictions)) (nreverse predictions))
      (cond
        ((and (notany #'(lambda(x) (alpha-char-p x)) line) (string/= "" aux))
         (push aux predictions) (setq aux "")) ;new-sequence
        ((position command-symbol line) nil) ;discard comments
        (t (setq aux (concatenate 'string aux (identify-region line))))))))


