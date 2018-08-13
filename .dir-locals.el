;;; Project-wide Emacs settings
;;; Directory Local Variables
;;; For more information see (info "(emacs) Directory Variables")

((nil
  (c-file-style . "bsd")
  (cmake-ide-cmake-opts . "-DBUILD_TESTS=ON -DBUILD_PYTHON_MODULE=ON")
  (cmake-ide-make-command . "make run_tests"))
 (c-mode
  (mode . c++)
  (flycheck-gcc-language-standard . "c++14")
  (flycheck-clang-language-standard . "c++14")
  (flycheck-clang-standard-library . "libc++")
  (indent-tabs-mode . nil)
  (c-basic-offset . 4))
 (c++-mode
  (flycheck-gcc-language-standard . "c++14")
  (flycheck-clang-language-standard . "c++14")
  (flycheck-clang-standard-library . "libc++")
  (indent-tabs-mode . nil)
  (c-basic-offset . 4)))
