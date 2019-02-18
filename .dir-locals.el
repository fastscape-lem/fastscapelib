;;; Project-wide Emacs settings
;;; Directory Local Variables
;;; For more information see (info "(emacs) Directory Variables")

((nil
  (c-file-style . "bsd")
  (cmake-ide-cmake-opts . "-DDOWNLOAD_GTEST=ON -DBUILD_TESTS=ON -DBUILD_PYTHON_MODULE=ON -DBUILD_BENCHMARK=ON"))
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
  (flycheck-clang-warnings . ("everything" "no-pragma-once-outside-header" "no-c++98-compat"))
  (indent-tabs-mode . nil)
  (c-basic-offset . 4)
  (c-file-offsets
   (innamespace . -))))
