%%
%%  A LaTeX2e class defining a UCL MRes, MPhil or PhD thesis.
%%
%%  This class was originally created (and is copyright) in 1996 by Russel Winder.
%%  It is maintained and distributed with permission, by Ian Kirker in
%%   the Research Computing team at ISD.
%%  (rc-support@ucl.ac.uk, or i.kirker@ucl.ac.uk)
%%
%%  This version differs from the commonly available version of thesis.cls, in
%%   being rebased on Russel's un-docstripped version from KCL, and thus includes,
%%   among other changes, significantly more documentation.
%%
%%  The definitive formatting reference is here:
%%   http://www.ucl.ac.uk/current-students/research_degrees/thesis_formatting 
%%
%%  Several quotes from that page are included at relevant points.
%%
%%  Many stylistic choices that used to be mandatory are now not specified, so
%%   there are many choices made here that no longer have to be set this way.
%%  Regardless, they are (in every case I've noticed) sound stylistic choices.
%%  It's a good idea to check the page above before changing something.

\def \fileversion{v1.6}
\def \filedate{2014/08/12}
\def \docdate{\filedate}

\newcommand \@doctor{Doctor of Philosophy}
\newcommand \@master{Master of Philosophy}
\newcommand \@mres{Master of Research}
\newcommand \@degree@string{\@doctor}

%  The usual structure of class preamble follows.
%
%  Assume A4 paper.  Ensure single column.  Assume double-sided as default.  Ensure A4 margins
%  and Times typeface (sans serif is Helvetica).  Assume 10pt as default.
%
%  Make the possibility of special activity in draft mode.  To do this create a boolean as a
%  flag and ensure that the appropriate options are passed onto the superclass.

\NeedsTeXFormat{LaTeX2e}[1997/12/01]
\ProvidesClass{thesis}[\filedate\space\fileversion\space UCL MRes/MPhil/PhD Thesis class]
\DeclareOption{mres}{\renewcommand\@degree@string{\@mres}}
\DeclareOption{MRes}{\renewcommand\@degree@string{\@mres}}
\DeclareOption{mphil}{\renewcommand\@degree@string{\@master}}
\DeclareOption{MPhil}{\renewcommand\@degree@string{\@master}}
\DeclareOption{phd}{\renewcommand\@degree@string{\@doctor}}
\DeclareOption{PhD}{\renewcommand\@degree@string{\@doctor}}
\DeclareOption{twocolumn}{\ClassError{thesis}{Option `twocolumn' not supported for UCL theses}{}}
\newif \if@draftmode
\@draftmodefalse
\DeclareOption{draft}{\@draftmodetrue\PassOptionsToClass{draft}{book}}
\DeclareOption{final}{\@draftmodefalse\PassOptionsToClass{final}{book}}
\PassOptionsToClass{a4paper}{book}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\ProcessOptions
\LoadClass{book}
\RequirePackage{mathptmx}

\setlength \hoffset{-25.4mm}
\setlength \voffset{-25.4mm}
\setlength \textwidth{140mm}
\setlength \textheight{232mm} % I think you need some extra vertical space for headers and page numbers
\setlength \topmargin{20mm}
\setlength \oddsidemargin {39mm} % Allow a mm for the bleed.
\if@twoside
  \setlength \evensidemargin {20mm}
\else
  \setlength \evensidemargin {39mm}
\fi

\def \ps@headings{%
    \def \@oddhead{\normalfont \rmfamily \slshape \hfill \rightmark \hfill \thepage}%
    \def \@evenhead{\normalfont \rmfamily \slshape \thepage \hfill \leftmark \hfill}%
    \let \@oddfoot = \@empty%
    \let \@evenfoot = \@empty%
    \let \@mkboth = \markboth%
    \def \chaptermark##1{%
        \markboth{%
            \ifnum \c@secnumdepth >\m@ne%
                \if@mainmatter%
                    \@chapapp\ \thechapter. \ %
                \fi%
            \fi%
            ##1%
        }{}%
    }%
    \def \sectionmark##1{%
        \markright {%
            \ifnum \c@secnumdepth >\z@%
                \thesection. \ %
            \fi%
            ##1%
        }%
    }%
}
\pagestyle{headings}
\renewcommand \footnoterule {%
    \kern -3\p@%
    \hrule%
    \kern 2.6\p@%
}

\setlength \parindent{2em}
\setlength \parskip  {\z@ \@plus 1\p@}
\setlength \headsep  {16\p@}
\ifcase \@ptsize
    \setlength \footskip{26\p@}
\or
    \setlength \footskip{27\p@}
\or
    \setlength \footskip{28\p@}
\fi
\setlength \partopsep {\z@}
\tolerance = \@m
\emergencystretch = 10\p@
\vfuzz = 2\p@
\hfuzz = 2\p@
\flushbottom
\renewcommand \baselinestretch{1.5}
\renewcommand \tableofcontents{%
    \if@twocolumn%
        \@restonecoltrue \onecolumn%
    \else%
        \@restonecolfalse%
    \fi%
    \chapter*{%
        \contentsname%
        \@mkboth{\contentsname}{\contentsname}%
    }%
    \@starttoc{toc}%
    \if@restonecol%
        \twocolumn%
    \fi%
}
\renewcommand \listoffigures{%
    \if@twocolumn%
        \@restonecoltrue \onecolumn%
    \else%
        \@restonecolfalse%
    \fi%
    \chapter*{%
        \listfigurename%
        \@mkboth{\listfigurename}{\listfigurename}%
    }%
    \@starttoc{lof}%
    \if@restonecol%
        \twocolumn%
    \fi%
}
\renewcommand \listoftables{%
    \if@twocolumn%
        \@restonecoltrue \onecolumn%
    \else%
        \@restonecolfalse%
    \fi%
    \chapter*{%
        \listtablename%
        \@mkboth{\listtablename}{\listtablename}%
    }%
    \@starttoc{lot}%
    \if@restonecol%
        \twocolumn%
    \fi%
}
\renewenvironment{thebibliography}[1]{%
    \chapter*{%
        \bibname%
        \@mkboth{\bibname}{\bibname}%
    }%

    \list{\@biblabel{\arabic{enumiv}}}%
         {\settowidth \labelwidth{\@biblabel{#1}}%
          \leftmargin \labelwidth%
          \advance \leftmargin \labelsep%
          \@openbib@code%
          \usecounter{enumiv}%
          \let \p@enumiv = \@empty%
          \renewcommand \theenumiv{\arabic{enumiv}}%
    }%
    \sloppy%
    \clubpenalty 4000%
    \@clubpenalty \clubpenalty%
    \widowpenalty 4000%
    \sfcode`\.=\@m%
}{%
    \def \@noitemerr{\@latex@warning{Empty `thebibliography' environment}}%
    \endlist%
}
\renewenvironment{theindex}{%
    \if@twocolumn%
        \@restonecolfalse%
    \else%
        \@restonecoltrue%
    \fi%
    \columnseprule \z@%
    \columnsep 35\p@%
    \twocolumn[\@makeschapterhead{\indexname}]%
    \@mkboth{\indexname}{\indexname}%

    \thispagestyle{plain}%
    \small%
    \parindent \z@%
    \parskip \z@ \@plus .3\p@%
    \relax%

    \let \item = \@idxitem%
}{%
    \if@restonecol%
        \onecolumn%
    \else%
        \clearpage%
    \fi%
}

\def\department#1{\gdef\@department{#1}}
\def\@department{\@latex@warning@no@line{No \noexpand\department given}}

\def\sponsor#1{\gdef\@sponsor{#1}}
\def\@sponsor{\@latex@warning@no@line{No \noexpand\sponsor given}}
\renewcommand \maketitle{%
    \setcounter{page}{1}%
    \thispagestyle{empty}%
    \@maketitle%
    \setcounter{footnote}{0}%
    \let \thanks = \relax%
    \gdef \@address{}%
    \gdef \@thanks{}%
    %\gdef \@author{}%
    \gdef \@department{}%
    \gdef \@sponsor{}%
    \gdef \@title{}%
    \let \maketitle = \relax%
}

\newcommand \@maketitle{%
    \newpage%
    \null%
    \vspace*{4em}%

    \begin{center}%
    \vspace{-2cm}
        {\huge \bfseries \@title}\\[5em]%
        {\Large \itshape \@author}\\%
    \end{center}%
    \vfill%
    \begin{center}%
    A dissertation submitted in fulfillment \\
    of the requirements for the degree of \\
    \vspace{0.5cm}
    \textbf{\@degree@string} \\
    at \\
    \textbf{King's College London}\\
    \@department
    \end{center}%
    \vspace{1em}%
    \begin{center}%
    \includegraphics[scale=0.125]{assets/kcl-crest.png}
    \end{center}%
    \vspace{1em}
    \begin{center}%
    funded by \textbf{\@sponsor}\\
    \vspace{1em}
    \includegraphics[scale=0.3]{assets/microsoft-logo.pdf}
    \end{center}%
    \vspace{1em}%
    \begin{center}%
    \@date%
    \end{center}%
    \if@twoside %
      \newpage%
      ~\\
      \newpage%
    \fi
}

       \newcommand \makedeclaration {%
         \clearpage%
         I, \@author , confirm that the work presented in this thesis is my own.
         Where information has been derived from other sources, I confirm that this has been indicated in the work.
         \gdef \@author{}%
         \clearpage%
       }%

\renewcommand \chapter{%
    \if@openright%
        \cleardoublepage%
    \else%
        \clearpage%
    \fi
    \thispagestyle{empty}%
    \global \@topnum \z@%
    \@afterindentfalse%
    \secdef \@chapter \@schapter%
}
\def \@makechapterhead#1{%
    \vspace*{40\p@}%
    {\parindent \z@%
     \raggedright%
     \reset@font%
     \ifnum \c@secnumdepth >\m@ne%
         \if@mainmatter%
             \Large \bfseries \@chapapp{} \thechapter%
             \par%
             \vskip 20\p@%
         \fi%
     \fi%
     \huge \bfseries #1%
     \par%
     \nobreak%
     \vskip 40\p@%
    }%
}
\def \@makeschapterhead#1{%
    \vspace*{50\p@}%
    {\parindent \z@%
     \raggedright%
     \reset@font%
     \huge \bfseries #1%
     \par%
     \nobreak%
     \vskip 40\p@%
    }%
}
\renewcommand \section{\@startsection {section}{1}{\z@}%
                                      {-1ex \@plus -.4ex \@minus -.2ex}%
                                      {.6ex \@plus .1ex}%
                                      {\reset@font \Large \bfseries}}
\renewcommand \subsection{\@startsection{subsection}{2}{\z@}%
                                        {-.7ex\@plus -.3ex \@minus -.2ex}%
                                        {.4ex \@plus .1ex}%
                                        {\reset@font \large \bfseries}}
\renewcommand \subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                        {-.5ex\@plus -.2ex \@minus -.2ex}%
                                        {.3ex \@plus .1ex}%
                                        {\reset@font \large}}
\renewcommand \paragraph{\@startsection{paragraph}{4}{\z@}%
                                       {.5ex \@plus .1ex \@minus .2ex}%
                                       {-.2em}%
                                       {\reset@font \normalsize \bfseries}}
\renewcommand \subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                          {-.3ex \@plus .1ex \@minus .2ex}%
                                          {-.2em}%
                                          {\reset@font \normalsize \bfseries}}
\setcounter{secnumdepth}{2}
\newcommand \@abstract{Abstract}
\newenvironment{abstract}{%
    \chapter*{\@abstract}%
    \@mkboth{\@abstract}{\@abstract}%
}{%
}

\newcommand \@acknowledgements{Acknowledgements}
\newenvironment{acknowledgements}{%
    \chapter*{\@acknowledgements}%
    \@mkboth{\@acknowledgements}{\@acknowledgements}%
}{%
}

\endinput