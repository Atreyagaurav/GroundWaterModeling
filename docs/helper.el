(defun wrap-in-src (beg end)
  (interactive "r")
  (save-excursion
    (goto-char end)
    (insert "#+END_SRC\n")
    (goto-char beg)
    (insert "#+BEGIN_SRC python\n")))

(global-set-key (kbd "<menu>") 'wrap-in-src)
