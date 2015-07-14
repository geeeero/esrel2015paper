\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{opacpstr}
              [2015/07/03 v0.1
               poster using TU/e 2008 corporate identity, based on mnposter]
% v1.1: updated itemize and enumerate in accordance with corporate identity
%
%
% let class by based on article.cls
%\LoadClass[12pt,twocolumn]{article}
\LoadClass[12pt]{article}
%page setup to A3, since that is the size for which the design of the
%``TU/e Wetenschapsposter'' has been defined.
\setlength\paperheight{420mm}%
\setlength\paperwidth{297mm}
\setlength{\hoffset}{0mm}
\setlength\oddsidemargin{13mm}
\addtolength\oddsidemargin{-1in}%to correct for standard 1in
\setlength\evensidemargin\oddsidemargin
\setlength\textwidth{273mm}
\setlength\columnsep{9mm}%NB columnwidth should be determined by LaTeX to 132mm
\setlength{\voffset}{0mm}
\setlength{\topmargin}{36.5mm}
\addtolength\topmargin{-1in}%to correct for standard 1in
\setlength{\headheight}{42mm}
\setlength{\headsep}{17mm}
\setlength\textheight{295mm}
\setlength\footskip{23mm}


%
% required before defining class options
\RequirePackage{color}
\ifx\pdfoutput\undefined
  \RequirePackage[dvips]{graphicx}
\else
  \RequirePackage[pdftex]{graphicx}
  \pdfoutput=1
  \pdfcompresslevel=9
  \pdfpagewidth=297mm
  \pdfpageheight=420mm
\fi
\RequirePackage{pgfpages}
\RequirePackage{url}


%
% class file options
\DeclareOption{red}{\def\themecolor{0}}
\DeclareOption{orange}{\def\themecolor{1}}
\DeclareOption{green}{\def\themecolor{2}}
\DeclareOption{cyan}{\def\themecolor{3}}
\DeclareOption{blue}{\def\themecolor{4}}
\DeclareOption{rood}{\def\themecolor{0}}
\DeclareOption{oranje}{\def\themecolor{1}}
\DeclareOption{groen}{\def\themecolor{2}}
\DeclareOption{cyaan}{\def\themecolor{3}}
\DeclareOption{middenblauw}{\def\themecolor{4}}
\DeclareOption{a4}{\pgfpagesuselayout{resize to}[a4paper]}
\DeclareOption{a1plus}{\pgfpagesphysicalpageoptions{logical pages=1,%
                                                    physical height=896.6mm,%
                                                    physical width=650.5mm}
                       \pgfpageslogicalpageoptions{1}{resized width=594mm,%
                                                      resized height=841mm,%
                                                      center=\pgfpoint{0.5\pgfphysicalwidth}{0.5\pgfphysicalheight},%
                                                      border code=\tueposterbordercode
                                                     }
                      }
\ExecuteOptions{red}
\ProcessOptions\relax

%
% (additionally) required packages
\RequirePackage[T1]{fontenc}
\RequirePackage{pgf}
%\RequirePackage[mtbold,mtpluscal]{mathtime}


% Define TU/e colors
% Main colors
\definecolor{tuepmsprocesscyan}{cmyk}{1,0,0,0}%PMS Process Cyan; RGB 0,162,222
\definecolor{tuecyan}{cmyk}{1,0,0,0}%PMS Process Cyan; RGB 0,162,222
\definecolor{tuepmswarmred}{cmyk}{0,1,0.9,0}%PMS Warm Red; RGB 247,49,49
\definecolor{tuepms206}{cmyk}{0,1,0.5,0}%PMS 206; RGB 214,0,74
\definecolor{tuepms226}{cmyk}{0.1,1,0,0}%PMS 226; RGB 214,0,123
\definecolor{tuepms253}{cmyk}{0.5,0.9,0,0}%PMS 253; RGB 173,32,173
\definecolor{tuepms300}{cmyk}{1,0.5,0,0}%PMS 300; RGB 0,102,204
\definecolor{tuepms2748}{cmyk}{1,0.8,0,0.2}%PMS 2748; RGB 16,16,115
% Support colors
\definecolor{tuepms137}{cmyk}{0,0.5,1,0}%PMS 137; RGB 255,154,0
\definecolor{tuepmsyellow012}{cmyk}{0,0.1,1,0}%PMS Yellow 012; RGB 255,221,0
\definecolor{tuepms396}{cmyk}{0.15,0,1,0}%PMS 396; RGB 206,223,0
\definecolor{tuepmsgreen}{cmyk}{1,0,0.6,0}%PMS Green; RGB 0,172,130
\definecolor{tuepms3135}{cmyk}{1,0.1,0.25,0}%PMS 3135; RGB 0,146,181
\definecolor{tuepms375}{cmyk}{0.55,0,1,0}%PMS 375; RGB 132,210,0
% different names for some colors
\definecolor{tuewarmred}{cmyk}{0,1,0.9,0}%PMS Warm Red; RGB 247,49,49
\definecolor{tuecyan}{cmyk}{1,0,0,0}%PMS Process Cyan; RGB 0,162,222
\definecolor{tuered}{cmyk}{0,1,0.5,0}%PMS 206; RGB 214,0,74
\definecolor{tueblue}{cmyk}{1,0.5,0,0}%PMS 300; RGB 0,102,204
\definecolor{tuedarkblue}{cmyk}{1,0.8,0,0.2}%PMS 2748; RGB 16,16,115
\definecolor{tueorange}{cmyk}{0,0.5,1,0}%PMS 137; RGB 255,154,0
\definecolor{tueyellow}{cmyk}{0,0.1,1,0}%PMS Yellow 012; RGB 255,221,0
\definecolor{tuelightgreen}{cmyk}{0.55,0,1,0}%PMS 375; RGB 132,210,0
\definecolor{tuegreen}{cmyk}{1,0,0.6,0}%PMS Green; RGB 0,172,130

%
% TU/e fonts
\renewcommand{\sfdefault}{phv}
\renewcommand{\rmdefault}{phv}
\renewcommand{\ttdefault}{pcr}
\renewcommand\normalsize{\@setfontsize\normalsize\@xiipt{17}}
\renewcommand\small{\@setfontsize\small{9.5}{14}}
\renewcommand\footnotesize{\@setfontsize\footnotesize\@viipt\@xiipt}
\renewcommand\scriptsize{\@setfontsize\scriptsize\@viiipt{9.5}}
\renewcommand\large{\@setfontsize\large{16}{17}}%used by section command...
\renewcommand\section{\@startsection {section}{1}{\z@}%
                                   {\baselineskip}
                                   {0.1pt}%
                                   {\color{tuecyan}\normalfont\large\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {0pt}%
                                     {-1ex}%
                                     {\normalfont\bfseries}}
% No numbering of sections, subsections, etc.
\setcounter{secnumdepth}{0}

%
%
% Define header (top of poster) and footer (bottom of poster)
% i.e., everything except for the two columns of text.
%
\newcommand{\ps@poster}{%
   \renewcommand{\@oddhead}{\rlap{\hskip-7mm
                                  \begin{pgfpicture}{0mm}{0mm}{285mm}{17mm}
                                     \pgftext[at=\pgfpoint{7mm}{54mm},left,bottom]{%
                                       \color{tuedarkblue}
                                       \normalsize\bfseries\@topline
                                     }
                                     \pgftext[at=\pgfpoint{285mm}{49mm},right,bottom]{%
                                       \includegraphics[width=107.5mm]{tuelogo}%
                                     }
                                     \pgfsetcolor{tuered}
                                     \pgfpathmoveto{\pgfpoint{0mm}{49mm}}
                                     \pgfpathlineto{\pgfpoint{285mm}{49mm}}
                                     \pgfusepath{fill}
                                     \ifcase\themecolor
                                       \pgfsetcolor{tuered}\or
                                       \pgfsetcolor{tueorange}\or
                                       \pgfsetcolor{tuelightgreen}\or
                                       \pgfsetcolor{tuecyan}\or
                                       \pgfsetcolor{tueblue}
                                     \fi
                                     \pgfpathmoveto{\pgfpoint{0mm}{-7mm}}
                                     \pgfpathlineto{\pgfpoint{188mm}{-7mm}}
                                     \pgfpathlineto{\pgfpoint{204mm}{49mm}}
                                     \pgfpathlineto{\pgfpoint{0mm}{49mm}}
                                     \pgfusepath{fill}
                                   \end{pgfpicture}}%
                            \begin{minipage}[b][42mm][t]{187mm}%
                              \rightskip=0pt plus2em
                              \parshape 4 0mm 187mm 0mm 184mm 0mm 180mm 0mm 176mm%
                              \looseness=-1%
                              \color{white}%
                              \fontsize{30pt}{36pt}\selectfont%
                              \@title\par
                              \vskip9pt
                              \ifcase\prevgraf
                                \or\parshape 5 0mm 185mm 0mm 183mm 0mm 181mm 0mm 179.5mm 0mm 176mm%regelbreedtes als titel 1 regel beslaat
                                \or\parshape 3 0mm 181mm 0mm 179.5mm 0mm 178mm%regelbreedtes als titel 2 regels beslaat
                                \or\parshape 2 0mm 177mm 0mm 176mm%regelbreedtes als titel 3 regels beslaat
                                \else\parshape 1 0mm 176mm%
                              \fi
                              \looseness=-1%
                              \normalsize
                              \@author\par
                            \end{minipage}\hfil
                            \llap{\begin{minipage}[b][42mm][t]{58mm}
                              \fontsize{8.5pt}{14pt}\selectfont
                              \raggedright
                              \@logo
                            \end{minipage}\hskip-5mm}%
                           }
   \renewcommand{\@evenhead}{\@oddhead}
   \renewcommand{\@oddfoot}{\rlap{\hskip-7mm
                                  \begin{pgfpicture}{0mm}{0mm}{285mm}{17mm}
                                     \pgfsetcolor{tuedarkblue}
                                     \pgfpathmoveto{\pgfpoint{0mm}{0mm}}
                                     \pgfpathlineto{\pgfpoint{285mm}{0mm}}
                                     \pgfpathlineto{\pgfpoint{285mm}{17mm}}
                                     \pgfpathlineto{\pgfpoint{0mm}{17mm}}
                                     \pgfusepath{fill}
                                   \end{pgfpicture}
                                  }
                            \begin{minipage}[b][17mm][c]{271mm}
                              \color{white}
                              \normalsize / \@department
                            \end{minipage}
                           }
   \renewcommand{\@evenfoot}{\@oddfoot}
}

%
% page setup
\pagestyle{poster}
\let\marginpar\undefined%undefine \marginpar
\setlength\parindent{0mm}
\setlength\parskip{0mm}
\setlength\leftmargini  {1.4em}
\setlength\leftmarginii \leftmargini
\setlength\leftmarginiii\leftmargini
\setlength\leftmarginiv \leftmargini
\setlength\leftmarginv  \leftmargini
\setlength\leftmarginvi \leftmargini
\setlength\labelsep\leftmargini
\renewcommand\labelitemi{\rlap{\color{tuecyan}\textbullet}}
\renewcommand\labelitemii{\rlap{\color{tuecyan}\normalfont\bfseries \textendash}}
\renewcommand\labelitemiii{\rlap{\color{tuecyan}\textasteriskcentered}}
\renewcommand\labelitemiv{\rlap{\color{tuecyan}\textperiodcentered}}
\renewcommand*\descriptionlabel[1]{\hspace\labelsep
                                \color{tuecyan}\normalfont\bfseries #1}
\def\enumerate{%
  \ifnum \@enumdepth >\thr@@\@toodeep\else
    \advance\@enumdepth\@ne
    \edef\@enumctr{enum\romannumeral\the\@enumdepth}%
      \expandafter
      \list
        \csname label\@enumctr\endcsname
        {\usecounter\@enumctr\def\makelabel##1{\hss\rlap{##1}}}%
  \fi}
\long\def\@makecaption#1#2{%
  \small%use small instead of normalsize
  \vskip\abovecaptionskip
  \sbox\@tempboxa{#1: #2}%
  \ifdim \wd\@tempboxa >\hsize
    #1: #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\box\@tempboxa\hfil}%left aligned instead of centered
  \fi
  \vskip\belowcaptionskip}
%
%set spacing of listings
\def\@listi{\leftmargin\leftmargini
            \labelwidth\leftmargini
            \advance\labelwidth-\labelsep
            \parsep \z@ \@plus2\p@ \@minus\p@
            \topsep \z@ \@plus2\p@ \@minus4\p@
            \itemsep\z@ \@plus2\p@ \@minus\p@
            \itemsep   \parsep}
\let\@listI\@listi
\@listi
\def\@listii {\leftmargin\leftmarginii
              \labelwidth\leftmarginii
              \advance\labelwidth-\labelsep
              \topsep    \z@ \@plus2\p@ \@minus\p@
              \parsep    \z@ \@plus\p@  \@minus\p@
              \itemsep   \parsep}
\def\@listiii{\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep
              \topsep    \z@ \@plus\p@\@minus\p@
              \parsep    \z@
              \partopsep \z@ \@plus\z@ \@minus\p@
              \itemsep   \parsep}



% define tick marks
\newcommand{\tueposterbordercode}{%
  \pgfusepath{discard}
  \pgfpathmoveto{\pgfpoint{0mm}{-2mm}}
  \pgfpathlineto{\pgfpoint{0mm}{-12mm}}
  \pgfstroke
  \pgfpathmoveto{\pgfpoint{-2mm}{0mm}}
  \pgfpathlineto{\pgfpoint{-12mm}{0mm}}
  \pgfstroke
  \pgfpathmoveto{\pgfpoint{0mm}{422mm}}
  \pgfpathlineto{\pgfpoint{0mm}{432mm}}
  \pgfstroke
  \pgfpathmoveto{\pgfpoint{-2mm}{420mm}}
  \pgfpathlineto{\pgfpoint{-12mm}{420mm}}
  \pgfstroke
  \pgfpathmoveto{\pgfpoint{297mm}{422mm}}
  \pgfpathlineto{\pgfpoint{297mm}{432mm}}
  \pgfstroke
  \pgfpathmoveto{\pgfpoint{299mm}{420mm}}
  \pgfpathlineto{\pgfpoint{309mm}{420mm}}
  \pgfstroke
  \pgfpathmoveto{\pgfpoint{297mm}{-2mm}}
  \pgfpathlineto{\pgfpoint{297mm}{-12mm}}
  \pgfstroke
  \pgfpathmoveto{\pgfpoint{299mm}{0mm}}
  \pgfpathlineto{\pgfpoint{309mm}{0mm}}
  \pgfstroke
}


\newcommand\witregel{\vskip\baselineskip}
\let\blankline\witregel


\def\title#1{\gdef\@title{#1}}
\title{}
\def\author#1{\gdef\@author{#1}}
\author{}
\def\department#1{\gdef\@department{#1}}
\department{school of industrial engineering}
\def\otherlogos#1{\gdef\@otherlogos{#1}}
\otherlogos{}

\def\topline#1{\gdef\@topline{#1}}
\def\logo#1{\gdef\@logo{#1}}
\topline{Operations, Planning, Accounting \& Control \quad \url{www.opac.ieis.tue.nl}}
\logo{%\hfil\includegraphics[height=20mm]{mnlogo}\hfil
      %\vskip2mm
      \vbox to 18mm{\vfil\hbox to 52mm{\hfil\@otherlogos\hfil}}%
     }


\raggedright


\endinput
