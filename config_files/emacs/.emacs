(set-language-environment "UTF-8")
(setq custom-file "~/.emacs-customize") ;; changes generated by the `customize` command are written to that file
(load "~/.emacs-customize")

(menu-bar-mode -1) ;; remove useless toolbars
(tool-bar-mode -1)
(scroll-bar-mode -1)
(column-number-mode t) ;; shows column numbers

(setq indent-tabs-mode nil) ;; spaces instead of tabs
(setq c-basic-offset 4) ;; 4 spaces for indent
(setq c-default-style "k&amp;r") ;; indentation style - may be "bsd", "whitesmith", "gnu" (default)

;; Shortcuts - see http://www.cs.rutgers.edu/LCSR-Computing/some-docs/emacs-chart.html
(global-set-key (kbd "C-c r") 'replace-string) ;; example
