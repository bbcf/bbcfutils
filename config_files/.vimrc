
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Various options                                             "
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" See the title in the terminal "
set title

" Lots of colors "
set t_Co=256

" No old vim compatibility "
set nocompatible

" Activate mouse "
set mouse=a

" See the line numbers "
set number

" See the col numbers "
set ruler

" Always see the status line "
set laststatus=2

" Highlight searches "
set hlsearch

" Highlight pairs of parenthesis "
set showmatch

" Highlight current line "
set cursorline

" Don't auto indent "
set noautoindent

" The tab key makes spaces "
set expandtab

" A tab is four spaces "
set tabstop=4

" An indent is four spaces "
set shiftwidth=4

" This is our standard "
set fileencoding=utf-8
set encoding=utf-8

" Show real tabs with icon "
set listchars=tab:⇥\
set list

" No special auto formating "
set formatoptions=

" Wrapping text "
set wrap

" Wrapping cursor "
set whichwrap=<,>,h,l,[,]

" Set no timeouts "
set timeoutlen=1

" Long history "
set history=1000

" Lots of tabs "
set tabpagemax=30


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Various other stuff                                         "
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Easier leader character "
let mapleader = ","
let maplocalleader = ","

" Color scheme for python "
let g:python_slow_sync = 1

" File type recognition "
filetype on
filetype plugin on

" Syntax highlighting "
syntax on


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Autocompletion                                              "
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
autocmd FileType python     set omnifunc=pythoncomplete#Complete
autocmd FileType javascript set omnifunc=javascriptcomplete#CompleteJS
autocmd FileType html       set omnifunc=htmlcomplete#CompleteTags
autocmd FileType css        set omnifunc=csscomplete#CompleteCSS
autocmd FileType xml        set omnifunc=xmlcomplete#CompleteTags
autocmd FileType php        set omnifunc=phpcomplete#CompletePHP
autocmd FileType c          set omnifunc=ccomplete#Complete


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Code folding                                                "
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
set foldmethod=indent   " fold based on indent
set foldnestmax=1       " deepest fold possible is 1 levels
set nofoldenable        " don't fold by default
set foldlevel=1         " this is just what I use


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Auto remove whitespace                                      "
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
function! <SID>StripTrailingWhitespaces()
    " prepare: save last search, and cursor position.
    let _s=@/
    let l = line(".")
    let c = col(".")
    " do the business:
    %s/\s\+$//e
    " clean up: restore previous search history, and cursor position
    let @/=_s
    call cursor(l, c)
endfunction
autocmd BufWritePre * :call <SID>StripTrailingWhitespaces()


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Syntax highlighting                                         "
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Change colors when changing modes "
colorscheme modecluen
autocmd InsertEnter * colorscheme modecluei
autocmd InsertLeave * colorscheme modecluen

" General "
highlight Comment    guifg=darkgray gui=NONE ctermfg=darkgray
highlight String     guifg=#ff0087           ctermfg=198
highlight Identifier guifg=black    gui=bold ctermfg=black cterm=bold
highlight Statement  guifg=darkblue          ctermfg=darkblue
highlight Number     guifg=#d75f00           ctermfg=166
highlight PreProc    guifg=#5f00ff           ctermfg=57


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Syntax matching                                             "
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" Schemes "
highlight StyleWarning guifg=#ffd75f guibg=#ffaf00 ctermfg=221 ctermbg=214
highlight StyleError                 guibg=#ff6600             ctermbg=214

" All tabs "
highlight def link AllTabs StyleError
autocmd Syntax * syntax match AllTabs /\t/ containedin=ALL
autocmd BufEnter * match StyleError /\t/

" Trailing spaces "
highlight def link ExtraWhitespace StyleWarning
autocmd Syntax * syntax match ExtraWhitespace /\s\+$\| \+\ze\t/ containedin=ALL
