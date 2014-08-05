#******************************************************************************
#                                                                             *
#              *    **   *  *   *                           *                 *
#             * *  *  *  *  *  * *                          *                 *
#            ***** *  *  *  * *****  **  ***  *  *  ** ***  ***               *
#            *   * *  *  *  * *   * *  * *  * *  * *   *  * *  *              *
#            *   * *  *  *  * *   * *  * *  * *  *   * *  * *  *              *
#            *   *  ** *  **  *   *  *** ***   *** **  ***  *  *              *
#                                      * *             *                      *
#                                    **  *             *                      *
#                                                                             *
#******************************************************************************
#                                                                             *
#  This file is part of AQUAgpusph, a free CFD program based on SPH.          *
#  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>               *
#                                                                             *
#  AQUAgpusph is free software: you can redistribute it and/or modify         *
#  it under the terms of the GNU General Public License as published by       *
#  the Free Software Foundation, either version 3 of the License, or          *
#  (at your option) any later version.                                        *
#                                                                             *
#  AQUAgpusph is distributed in the hope that it will be useful,              *
#  but WITHOUT ANY WARRANTY; without even the implied warranty of             *
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
#  GNU General Public License for more details.                               *
#                                                                             *
#  You should have received a copy of the GNU General Public License          *
#  along with AQUAgpusph.  If not, see <http://www.gnu.org/licenses/>.        *
#                                                                             *
#******************************************************************************

import os


PURGE = ['.git', '.gitignore']
ICONS = {'folder' : r'\includegraphics[width=0.05\textwidth]{tree/Actions-document-open-folder-icon}',
         'file' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-x-office-document-icon}',
         'Releases notes' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-plain-icon}',
         'LICENSE' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-plain-icon}',
         'Code styling' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-plain-icon}',
         '.md' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-plain-icon}',
         'TODO' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-plain-icon}',
         '.cmake' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-plain-icon}',
         '.txt' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-plain-icon}',
         '.in' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-plain-icon}',
         '.log' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-plain-icon}',
         '.dat' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-plain-icon}',
         '.pdf' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-application-pdf-icon}',
         '.tex' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-x-bibtex-icon}',
         '.png' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-application-x-egon-icon}',
         '.jpg' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-application-x-egon-icon}',
         '.svg' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-application-x-egon-icon}',
         '.py' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-x-python-icon}',
         '.c' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-x-csrc-icon}',
         '.cl' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-x-csrc-icon}',
         '.cpp' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-x-c-plus-plus-src-icon}',
         '.h' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-x-chdr-icon}',
         '.hcl' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-text-x-chdr-icon}',
         '.xml' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-application-xml-icon}',
         '.sh' : r'\includegraphics[width=0.05\textwidth]{tree/Mimetypes-application-x-shellscript-icon}',
}


def list_files(startpath):
    txt = '\dirtree{%\n'
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        name = os.path.basename(root)
        purged = False
        for p in PURGE:
            if p in root:
                purged = True
                break
        if purged:
            continue
        txt += '.{} {{ {} {} }}.\n'.format(level + 1,
                                           ICONS['folder'],
                                           name)
        for f in files:
            if f in PURGE:
                continue
            icon = ICONS['file']
            for extension in ICONS.keys()[2:]:
                if f.endswith(extension):
                    icon = ICONS[extension]
                    break
            txt += '.{} {{ {} {} }}.\n'.format(level + 2,
                                               icon,
                                               f)
    txt += '}\n'
    # Fix some special characters
    txt = txt.replace('_', '\_')
    return txt

if __name__ == '__main__':
    script_path = os.path.dirname(os.path.realpath(__file__))
    root_path = os.path.abspath(os.path.join(script_path, '../../../../'))
    print(list_files(root_path))

