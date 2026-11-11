# -*- coding: utf-8 -*-
"""Build the guide to the scripts and the workflow (.docx), figures included.

A DERIVED document, like the cookbook: edit this script, re-run it, and
commit both. Editing the .docx by hand means the next regeneration discards
the change.

    python docs/make_scripts_guide.py                # beside this script
    python docs/make_scripts_guide.py --figures DIR  # the figures as PNG too

The .docx is written with the standard library alone (a .docx is a ZIP of
XML): python-docx, which docs/make_cookbook.py uses, is not in the model
environment of the server (2026-11-11), and matplotlib -- already a
dependency of the model -- draws the three workflow figures.
Set MM_GUIDE_OUT to write the document elsewhere.
"""
import io
import os
import re
import sys
import zipfile
from xml.sax.saxutils import escape as _esc

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get('MM_GUIDE_OUT',
                     os.path.join(HERE, 'MARMITES_scripts_guide.docx'))
DATE = '2026-11-11'
TITLE = 'MARMITES-MODFLOW 6: guide to the scripts and the workflow'

TEXT_W = 9638                  # dxa: A4 less 2 cm margins (17 cm)
FIG_W = 16.8                   # cm, a figure's width on the page
ACCENT = '1F4E79'
GREY = '595959'

# ===================================================================== #
#  a .docx with the standard library
# ===================================================================== #

NS = ('xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
      'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/'
      'relationships" '
      'xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/'
      'wordprocessingDrawing" '
      'xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
      'xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture"')
W_NS = 'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"'
XML = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
_INLINE = re.compile(r'(\*\*.+?\*\*|`.+?`)')


def _run(text, rpr=''):
    return '<w:r>%s<w:t xml:space="preserve">%s</w:t></w:r>' % (
        '<w:rPr>%s</w:rPr>' % rpr if rpr else '', _esc(text))


def runs(text):
    """Inline markup to runs: **bold** and `code`, a `code` inside a
    **bold** being both."""
    out = []
    for part in _INLINE.split(text):
        if not part:
            continue
        if part.startswith('**'):
            for k, piece in enumerate(part[2:-2].split('`')):
                if piece:       # odd pieces sit between two backticks
                    out.append(_run(piece, ('<w:rStyle w:val="CodeChar"/>'
                                            if k % 2 else '') + '<w:b/>'))
        elif part.startswith('`'):
            out.append(_run(part[1:-1], '<w:rStyle w:val="CodeChar"/>'))
        else:
            out.append(_run(part))
    return ''.join(out)


def para(text, style=None, ppr=''):
    """A paragraph; ``ppr``: paragraph properties after pStyle, in the
    schema's order (keepNext, numPr, spacing, ind, jc)."""
    return '<w:p><w:pPr>%s%s</w:pPr>%s</w:p>' % (
        '<w:pStyle w:val="%s"/>' % style if style else '', ppr, runs(text))


class Doc:
    def __init__(self):
        self.body = []
        self.media = []              # (name, bytes)
        self.nums = []               # abstractNumId of each w:num, numId 1..
        self._bm = 0

    # ---- text ---------------------------------------------------------
    def h(self, level, text, anchor=None):
        inner = runs(text)
        if anchor:
            self._bm += 1
            inner = ('<w:bookmarkStart w:id="%d" w:name="%s"/>%s'
                     '<w:bookmarkEnd w:id="%d"/>'
                     % (self._bm, anchor, inner, self._bm))
        self.body.append('<w:p><w:pPr><w:pStyle w:val="Heading%d"/></w:pPr>%s'
                         '</w:p>' % (level, inner))

    def p(self, text, style=None):
        self.body.append(para(text, style))

    def note(self, text):
        self.body.append(para(text, 'Note'))

    def code(self, lines):
        for k, ln in enumerate(lines):
            self.body.append(
                '<w:p><w:pPr><w:pStyle w:val="CodeBlock"/>%s</w:pPr>%s</w:p>'
                % ('<w:keepNext/>' if k < len(lines) - 1 else '',
                   _run(ln) if ln else ''))
        self.body.append(para('', 'Spacer'))

    def _list(self, items, abstract):
        self.nums.append(abstract)
        num = len(self.nums)
        for it in items:
            lvl, text = (it if isinstance(it, tuple) else (0, it))
            self.body.append(para(
                text, 'ListParagraph',
                '<w:numPr><w:ilvl w:val="%d"/><w:numId w:val="%d"/></w:numPr>'
                % (lvl, num)))

    def bullets(self, items):
        self._list(items, 0)

    def steps(self, items):
        self._list(items, 1)

    def link_list(self, items):
        """Contents: (anchor, text) lines linking to the headings."""
        for anchor, text in items:
            self.body.append(
                '<w:p><w:pPr><w:pStyle w:val="Contents"/></w:pPr>'
                '<w:hyperlink w:anchor="%s" w:history="1">%s</w:hyperlink>'
                '</w:p>' % (anchor, _run(text, '<w:rStyle w:val="Hyperlink"/>')))

    def pagebreak(self):
        self.body.append('<w:p><w:r><w:br w:type="page"/></w:r></w:p>')

    # ---- tables -------------------------------------------------------
    def table(self, head, rows, widths):
        """``widths``: relative column widths; cells take the inline
        markup, and a list of strings gives a cell several paragraphs."""
        tot = float(sum(widths))
        cols = [int(round(TEXT_W * w / tot)) for w in widths]
        cols[-1] += TEXT_W - sum(cols)
        border = ''.join('<w:%s w:val="single" w:sz="4" w:space="0" '
                         'w:color="BFBFBF"/>' % s for s in
                         ('top', 'left', 'bottom', 'right', 'insideH',
                          'insideV'))
        x = ['<w:tbl><w:tblPr><w:tblW w:w="%d" w:type="dxa"/><w:tblBorders>%s'
             '</w:tblBorders><w:tblLayout w:type="fixed"/><w:tblCellMar>'
             '<w:top w:w="50" w:type="dxa"/><w:left w:w="90" w:type="dxa"/>'
             '<w:bottom w:w="50" w:type="dxa"/><w:right w:w="90" '
             'w:type="dxa"/></w:tblCellMar></w:tblPr><w:tblGrid>%s'
             '</w:tblGrid>' % (TEXT_W, border, ''.join(
                 '<w:gridCol w:w="%d"/>' % c for c in cols))]

        def row(cells, header):
            out = ['<w:tr><w:trPr><w:cantSplit/>%s</w:trPr>'
                   % ('<w:tblHeader/>' if header else '')]
            for c, w in zip(cells, cols):
                paras = c if isinstance(c, (list, tuple)) else [c]
                out.append(
                    '<w:tc><w:tcPr><w:tcW w:w="%d" w:type="dxa"/>%s</w:tcPr>'
                    '%s</w:tc>' % (
                        w, '<w:shd w:val="clear" w:color="auto" '
                        'w:fill="DCE6F1"/>' if header else '',
                        ''.join(para(t, 'TableHead' if header else
                                     'TableText') for t in paras)))
            out.append('</w:tr>')
            return ''.join(out)
        x.append(row(head, True))
        x += [row(r, False) for r in rows]
        x.append('</w:tbl>')
        self.body.append(''.join(x))
        self.body.append(para('', 'Spacer'))

    # ---- figures ------------------------------------------------------
    def figure(self, png, w_cm, h_cm, caption):
        self.media.append(('image%d.png' % (len(self.media) + 1), png))
        n = len(self.media)
        cx, cy = int(w_cm * 360000), int(h_cm * 360000)
        self.body.append(
            '<w:p><w:pPr><w:pStyle w:val="Figure"/></w:pPr><w:r><w:drawing>'
            '<wp:inline distT="0" distB="0" distL="0" distR="0">'
            '<wp:extent cx="%d" cy="%d"/><wp:docPr id="%d" name="Figure %d"/>'
            '<wp:cNvGraphicFramePr><a:graphicFrameLocks noChangeAspect="1"/>'
            '</wp:cNvGraphicFramePr><a:graphic><a:graphicData uri="http://'
            'schemas.openxmlformats.org/drawingml/2006/picture"><pic:pic>'
            '<pic:nvPicPr><pic:cNvPr id="%d" name="image%d.png"/><pic:cNvPicPr/>'
            '</pic:nvPicPr><pic:blipFill><a:blip r:embed="rIdImg%d"/>'
            '<a:stretch><a:fillRect/></a:stretch></pic:blipFill><pic:spPr>'
            '<a:xfrm><a:off x="0" y="0"/><a:ext cx="%d" cy="%d"/></a:xfrm>'
            '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr>'
            '</pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing>'
            '</w:r></w:p>' % (cx, cy, n, n, n, n, n, cx, cy))
        self.body.append(para(caption, 'Caption'))

    # ---- the package --------------------------------------------------
    def _styles(self):
        def pstyle(sid, name, ppr='', rpr='', based='Normal', nxt=None,
                   q=True):
            return ('<w:style w:type="paragraph" w:styleId="%s"><w:name '
                    'w:val="%s"/><w:basedOn w:val="%s"/>%s%s%s%s</w:style>'
                    % (sid, name, based,
                       '<w:next w:val="%s"/>' % nxt if nxt else '',
                       '<w:qFormat/>' if q else '',
                       '<w:pPr>%s</w:pPr>' % ppr if ppr else '',
                       '<w:rPr>%s</w:rPr>' % rpr if rpr else ''))

        def size(pt):
            return '<w:sz w:val="%d"/><w:szCs w:val="%d"/>' % (2 * pt, 2 * pt)
        col = '<w:color w:val="%s"/>' % ACCENT
        s = [XML, '<w:styles %s><w:docDefaults><w:rPrDefault><w:rPr>'
             '<w:rFonts w:ascii="Calibri" w:hAnsi="Calibri" '
             'w:eastAsia="Calibri" w:cs="Calibri"/><w:sz w:val="21"/>'
             '<w:szCs w:val="21"/><w:lang w:val="en-GB"/></w:rPr>'
             '</w:rPrDefault><w:pPrDefault><w:pPr><w:spacing w:after="120" '
             'w:line="264" w:lineRule="auto"/></w:pPr></w:pPrDefault>'
             '</w:docDefaults>' % W_NS,
             '<w:style w:type="paragraph" w:default="1" w:styleId="Normal">'
             '<w:name w:val="Normal"/><w:qFormat/></w:style>',
             '<w:style w:type="character" w:default="1" '
             'w:styleId="DefaultParagraphFont"><w:name w:val="Default '
             'Paragraph Font"/><w:uiPriority w:val="1"/><w:semiHidden/>'
             '</w:style>',
             '<w:style w:type="table" w:default="1" w:styleId="TableNormal">'
             '<w:name w:val="Normal Table"/><w:uiPriority w:val="99"/>'
             '<w:semiHidden/><w:tblPr><w:tblInd w:w="0" w:type="dxa"/>'
             '<w:tblCellMar><w:top w:w="0" w:type="dxa"/><w:left w:w="108" '
             'w:type="dxa"/><w:bottom w:w="0" w:type="dxa"/><w:right '
             'w:w="108" w:type="dxa"/></w:tblCellMar></w:tblPr></w:style>',
             pstyle('Title', 'Title', '<w:spacing w:before="0" w:after="60"/>',
                    '<w:b/>' + col + size(24), nxt='Normal'),
             pstyle('Subtitle', 'Subtitle', '<w:spacing w:after="60"/>',
                    '<w:color w:val="%s"/>' % GREY + size(13), nxt='Normal'),
             pstyle('Heading1', 'heading 1', '<w:keepNext/><w:keepLines/>'
                    '<w:spacing w:before="360" w:after="120"/><w:outlineLvl '
                    'w:val="0"/>', '<w:b/>' + col + size(16), nxt='Normal'),
             pstyle('Heading2', 'heading 2', '<w:keepNext/><w:keepLines/>'
                    '<w:spacing w:before="260" w:after="100"/><w:outlineLvl '
                    'w:val="1"/>', '<w:b/>' + col + size(13), nxt='Normal'),
             pstyle('Heading3', 'heading 3', '<w:keepNext/><w:keepLines/>'
                    '<w:spacing w:before="200" w:after="80"/><w:outlineLvl '
                    'w:val="2"/>', '<w:b/>' + col + size(11), nxt='Normal'),
             pstyle('ListParagraph', 'List Paragraph',
                    '<w:spacing w:after="70"/>'),
             pstyle('Contents', 'Contents',
                    '<w:spacing w:after="40"/><w:ind w:left="284"/>', q=False),
             pstyle('TableText', 'Table Text', '<w:spacing w:after="30" '
                    'w:line="240" w:lineRule="auto"/>', size(9), q=False),
             pstyle('TableHead', 'Table Head', '<w:keepNext/><w:spacing '
                    'w:after="30" w:line="240" w:lineRule="auto"/>',
                    '<w:b/>' + size(9), q=False),
             pstyle('CodeBlock', 'Code Block', '<w:keepLines/><w:shd '
                    'w:val="clear" w:color="auto" w:fill="F2F2F2"/><w:spacing '
                    'w:after="0" w:line="240" w:lineRule="auto"/><w:ind '
                    'w:left="170" w:right="170"/>',
                    '<w:rFonts w:ascii="Consolas" w:hAnsi="Consolas" '
                    'w:cs="Consolas"/><w:sz w:val="17"/><w:szCs w:val="17"/>',
                    q=False),
             pstyle('Note', 'Note', '<w:pBdr><w:left w:val="single" w:sz="18" '
                    'w:space="8" w:color="%s"/></w:pBdr><w:shd w:val="clear" '
                    'w:color="auto" w:fill="EEF3F8"/><w:spacing w:before="80" '
                    'w:after="160"/><w:ind w:left="200" w:right="120"/>'
                    % ACCENT, q=False),
             pstyle('Figure', 'Figure', '<w:keepNext/><w:spacing '
                    'w:before="120" w:after="60"/><w:jc w:val="center"/>',
                    q=False),
             pstyle('Caption', 'caption', '<w:spacing w:after="200"/><w:jc '
                    'w:val="center"/>', '<w:i/><w:color w:val="%s"/>' % GREY
                    + size(9)),
             pstyle('Spacer', 'Spacer', '<w:spacing w:after="0" w:line="120" '
                    'w:lineRule="auto"/>', size(4), q=False),
             pstyle('Footer', 'footer', '<w:spacing w:after="0"/><w:jc '
                    'w:val="center"/>', '<w:color w:val="%s"/>' % GREY
                    + size(8), q=False),
             '<w:style w:type="character" w:styleId="CodeChar"><w:name '
             'w:val="Code Char"/><w:rPr><w:rFonts w:ascii="Consolas" '
             'w:hAnsi="Consolas" w:cs="Consolas"/><w:sz w:val="18"/>'
             '<w:szCs w:val="18"/></w:rPr></w:style>',
             '<w:style w:type="character" w:styleId="Hyperlink"><w:name '
             'w:val="Hyperlink"/><w:rPr><w:color w:val="%s"/><w:u '
             'w:val="single"/></w:rPr></w:style>' % ACCENT,
             '</w:styles>']
        return ''.join(s)

    def _numbering(self):
        def lvl(i, fmt, text, left):
            return ('<w:lvl w:ilvl="%d"><w:start w:val="1"/><w:numFmt '
                    'w:val="%s"/><w:lvlText w:val="%s"/><w:lvlJc w:val="left"/>'
                    '<w:pPr><w:ind w:left="%d" w:hanging="260"/></w:pPr>'
                    '</w:lvl>' % (i, fmt, text, left))
        s = [XML, '<w:numbering %s>' % W_NS,
             '<w:abstractNum w:abstractNumId="0"><w:multiLevelType '
             'w:val="hybridMultilevel"/>%s%s</w:abstractNum>'
             % (lvl(0, 'bullet', '•', 500), lvl(1, 'bullet', '–',
                                                      900)),
             '<w:abstractNum w:abstractNumId="1"><w:multiLevelType '
             'w:val="hybridMultilevel"/>%s%s</w:abstractNum>'
             % (lvl(0, 'decimal', '%1.', 500), lvl(1, 'bullet', '–',
                                                   900))]
        for k, a in enumerate(self.nums):
            s.append('<w:num w:numId="%d"><w:abstractNumId w:val="%d"/>'
                     '<w:lvlOverride w:ilvl="0"><w:startOverride w:val="1"/>'
                     '</w:lvlOverride></w:num>' % (k + 1, a))
        s.append('</w:numbering>')
        return ''.join(s)

    def save(self, path):
        doc = (XML + '<w:document %s><w:body>%s<w:sectPr><w:footerReference '
               'w:type="default" r:id="rIdFooter"/><w:pgSz w:w="11906" '
               'w:h="16838"/><w:pgMar w:top="1134" w:right="1134" '
               'w:bottom="1134" w:left="1134" w:header="567" w:footer="567" '
               'w:gutter="0"/></w:sectPr></w:body></w:document>'
               % (NS, ''.join(self.body)))
        footer = (XML + '<w:ftr %s><w:p><w:pPr><w:pStyle w:val="Footer"/>'
                  '</w:pPr>%s<w:r><w:fldChar w:fldCharType="begin"/></w:r>'
                  '<w:r><w:instrText xml:space="preserve"> PAGE </w:instrText>'
                  '</w:r><w:r><w:fldChar w:fldCharType="separate"/></w:r>'
                  '<w:r><w:t>1</w:t></w:r><w:r><w:fldChar w:fldCharType="end"/>'
                  '</w:r></w:p></w:ftr>'
                  % (W_NS, _run('%s  ·  %s  ·  page ' % (TITLE, DATE))))
        rel = ('<Relationship Id="%s" Type="http://schemas.openxmlformats.org/'
               'officeDocument/2006/relationships/%s" Target="%s"/>')
        doc_rels = [rel % ('rIdStyles', 'styles', 'styles.xml'),
                    rel % ('rIdNum', 'numbering', 'numbering.xml'),
                    rel % ('rIdSettings', 'settings', 'settings.xml'),
                    rel % ('rIdFooter', 'footer', 'footer1.xml')]
        doc_rels += [rel % ('rIdImg%d' % (k + 1), 'image', 'media/' + n)
                     for k, (n, _b) in enumerate(self.media)]
        rels_ns = ('xmlns="http://schemas.openxmlformats.org/package/2006/'
                   'relationships"')
        ct = 'application/vnd.openxmlformats-officedocument.wordprocessingml.'
        parts = [
            ('[Content_Types].xml', XML + '<Types xmlns="http://schemas.'
             'openxmlformats.org/package/2006/content-types"><Default '
             'Extension="rels" ContentType="application/vnd.openxmlformats-'
             'package.relationships+xml"/><Default Extension="xml" '
             'ContentType="application/xml"/><Default Extension="png" '
             'ContentType="image/png"/>' + ''.join(
                 '<Override PartName="/word/%s.xml" ContentType="%s%s+xml"/>'
                 % (n, ct, t) for n, t in (
                     ('document', 'document.main'), ('styles', 'styles'),
                     ('numbering', 'numbering'), ('settings', 'settings'),
                     ('footer1', 'footer'))) +
             '<Override PartName="/docProps/core.xml" ContentType="application'
             '/vnd.openxmlformats-package.core-properties+xml"/></Types>'),
            ('_rels/.rels', XML + '<Relationships %s>%s<Relationship '
             'Id="rIdCore" Type="http://schemas.openxmlformats.org/package/'
             '2006/relationships/metadata/core-properties" Target="docProps/'
             'core.xml"/></Relationships>'
             % (rels_ns, rel % ('rIdDoc', 'officeDocument',
                                'word/document.xml'))),
            ('docProps/core.xml', XML + '<cp:coreProperties xmlns:cp="http://'
             'schemas.openxmlformats.org/package/2006/metadata/core-'
             'properties" xmlns:dc="http://purl.org/dc/elements/1.1/" '
             'xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://'
             'www.w3.org/2001/XMLSchema-instance"><dc:title>%s</dc:title>'
             '<dc:creator>MARMITES project</dc:creator><dcterms:created '
             'xsi:type="dcterms:W3CDTF">%sT00:00:00Z</dcterms:created>'
             '</cp:coreProperties>' % (_esc(TITLE), DATE)),
            ('word/document.xml', doc),
            ('word/_rels/document.xml.rels', XML + '<Relationships %s>%s'
             '</Relationships>' % (rels_ns, ''.join(doc_rels))),
            ('word/styles.xml', self._styles()),
            ('word/numbering.xml', self._numbering()),
            ('word/settings.xml', XML + '<w:settings %s><w:zoom '
             'w:percent="100"/><w:defaultTabStop w:val="720"/><w:compat>'
             '<w:compatSetting w:name="compatibilityMode" w:uri="http://'
             'schemas.microsoft.com/office/word" w:val="15"/></w:compat>'
             '</w:settings>' % W_NS),
            ('word/footer1.xml', footer)]
        parts += [('word/media/' + n, b) for n, b in self.media]
        with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
            for name, data in parts:
                # a fixed date: the same content gives the same file
                zi = zipfile.ZipInfo(name, (2026, 11, 11, 0, 0, 0))
                zi.compress_type = zipfile.ZIP_DEFLATED
                z.writestr(zi, data if isinstance(data, bytes)
                           else data.encode('utf-8'))


# ===================================================================== #
#  the figures
# ===================================================================== #

INK = '#1A1A1A'
KIND = {                       # fill, edge
    'you': ('#DCE9F5', '#1F4E79'),      # what the modeller does (an app page)
    'script': ('#FFFFFF', '#1F4E79'),   # a script
    'data': ('#F0F0F0', '#8A8A8A'),     # files and folders
    'pest': ('#FBEBDD', '#B5651D'),     # a step of the calibration
    'frame': ('#FFFFFF', '#B5651D'),    # a container of steps
    'frame2': ('#FFFFFF', '#1F4E79'),
}
MONO = 'DejaVu Sans Mono'


def _canvas(w, h):
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(w / 2.54, h / 2.54))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w)
    ax.set_ylim(0, h)
    ax.axis('off')
    return fig, ax


def _box(ax, x, y, w, h, kind, title=None, lines=(), fs=6.6, pad=0.16,
         dashed=False):
    """A rounded box, its bold title and its lines; a line starting with
    '`' is a script name (monospace), and what follows a '|' in it is its
    description, in the text face on the same line."""
    from matplotlib.patches import FancyBboxPatch
    fc, ec = KIND[kind]
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle='round,pad=0,rounding_size=0.11', fc=fc,
        ec=ec, lw=0.9, ls='--' if dashed else '-', zorder=2))
    lh = fs * 0.0353 * 1.5
    ty = y + h - pad
    if title:
        ax.text(x + pad, ty, title, fontsize=fs + 0.6, fontweight='bold',
                color=INK, va='top', ha='left', zorder=3)
        ty -= lh * 1.12
    for ln in lines:
        mono = ln.startswith('`')
        name, _sep, what = (ln[1:] if mono else ln).partition('|')
        ax.text(x + pad, ty, name, fontsize=fs - (0.5 if mono else 0.0),
                color=INK, va='top', ha='left',
                family=MONO if mono else 'DejaVu Sans', zorder=3)
        if what:
            # a monospace character is 0.602 em wide
            dx = len(name) * (fs - 0.5) * 0.0353 * 0.602 + 0.18
            ax.text(x + pad + dx, ty, what, fontsize=fs, color=INK, va='top',
                    ha='left', zorder=3)
        ty -= lh


def _arrow(ax, p0, p1, rad=0.0, color='#404040', lw=0.9, style='-|>'):
    ax.annotate('', xy=p1, xytext=p0, zorder=4, arrowprops=dict(
        arrowstyle=style, lw=lw, color=color, shrinkA=0, shrinkB=0,
        mutation_scale=8, connectionstyle='arc3,rad=%g' % rad))


def _legend(ax, x, y, kinds, fs=6.4):
    from matplotlib.patches import FancyBboxPatch
    for kind, label in kinds:
        fc, ec = KIND[kind]
        ax.add_patch(FancyBboxPatch((x, y), 0.42, 0.26, fc=fc, ec=ec, lw=0.8,
                                    boxstyle='round,pad=0,rounding_size=0.05'))
        ax.text(x + 0.55, y + 0.13, label, fontsize=fs, color=INK,
                va='center', ha='left')
        x += 0.75 + 0.118 * len(label)


def _png(fig):
    import matplotlib.pyplot as plt
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=220)
    plt.close(fig)
    return buf.getvalue()


def fig_workflow():
    """Figure 1: from the data to a calibrated model -- what the modeller
    does, the script behind it, what it reads or writes."""
    W, H = 17.0, 18.1
    fig, ax = _canvas(W, H)
    x1, w1, x2, w2, x3, w3 = 0.95, 4.75, 6.25, 5.55, 12.35, 4.45
    rh, gap = 1.5, 0.34
    top = H - 0.95
    for x, w, t in ((x1, w1, 'What you do (app page)'),
                    (x2, w2, 'The script behind it'),
                    (x3, w3, 'What it reads or writes')):
        ax.text(x + w / 2, H - 0.42, t, fontsize=7.6, fontweight='bold',
                color=INK, ha='center', va='center')
    rows = [
        ('you', 'Start the app', [],
         ['`tools/launch_app.py', '`app/Home.py'],
         ['a local copy of code/,', 'when the repository is on',
          'a network drive']),
        ('you', 'Home', ['where things are,', 'on this machine'],
         ['`mm_paths.py'],
         ['`configs/paths.local.toml']),
        ('you', '1  Grid', ['GIS layers, conversion,', 'mesh'],
         ['`tools/gis_to_dataset.py', '`marmites_meshes.py'],
         ['GIS folder  →', '`example/<case>/', 'the mesh, cached']),
        ('you', '2 to 6', ['surface, soil, groundwater,',
                           'observations, plots'],
         ['`app/lib/schema.py, editor.py', '`marmites_config.py'],
         ['`configs/<model>.toml']),
        ('you', '7  Run', ['validate, launch,', 'follow the log'],
         ['`app/lib/checks.py, runs.py', '`tests/run_lamata_mf6.py',
          'the run driver (Figure 2)'],
         ['`runs/<id>/', '`MF6_ws_<grid>/', '`out_<stamp>_<tag>/']),
        ('you', '10  Results', ['figures, maps of one day,',
                                'satellite checks'],
         ['`app/lib/results.py', '`ppMF6/marmites_daymap.py',
          '`marmites_wetimage.py, _smrs.py'],
         ['`out_<stamp>_<tag>/_input', '`out_<stamp>_<tag>/_output']),
        ('pest', '8  PEST++ set-up', ['define, set up, launch'],
         ['`tools/pest_lamata_setup.py', '`tools/pest_lamata_run.py',
          'the calibration (Figure 3)'],
         ['PEST workspace:', '`org/  template/', '`master/  workers/']),
        ('pest', '9  PEST++ progress', ['follow the iterations'],
         ['`app/lib/pest_progress.py', '`tools/pest_lamata_report.py'],
         ['the master\'s record files', '`master/report_<k>/']),
        ('pest', '8  Export, then 7  Run', ['the calibrated model'],
         ['`tools/pest_lamata_export.py'],
         ['`configs/<model>_PEST.toml', 'then a run, as any other']),
    ]
    ys = []
    for k, (kind, title, sub, scripts, data) in enumerate(rows):
        y = top - (k + 1) * rh - k * gap
        ys.append(y)
        _box(ax, x1, y, w1, rh, kind, title, sub)
        _box(ax, x2, y, w2, rh, 'script', None, scripts)
        _box(ax, x3, y, w3, rh, 'data', None, data)
        ym = y + rh / 2
        _arrow(ax, (x1 + w1, ym), (x2, ym), color='#7A7A7A', lw=0.7)
        _arrow(ax, (x2 + w2, ym), (x3, ym), color='#7A7A7A', lw=0.7)
        if k:
            _arrow(ax, (x1 + w1 / 2, y + rh + gap), (x1 + w1 / 2, y + rh))
    # the calibrated model goes back to the Run page
    xa = 0.42
    ax.plot([x1, xa, xa], [ys[8] + rh / 2, ys[8] + rh / 2, ys[4] + rh / 2],
            color='#B5651D', lw=0.9, zorder=1)
    _arrow(ax, (xa, ys[4] + rh / 2), (x1, ys[4] + rh / 2), color='#B5651D')
    ax.text(xa - 0.12, (ys[8] + ys[4]) / 2 + rh / 2, 'run the calibrated model',
            fontsize=6.2, color=INK, rotation=90, va='center', ha='center')
    _legend(ax, x1, 0.22, [('you', 'a step on an app page'),
                           ('pest', 'a step of the calibration'),
                           ('script', 'scripts (paths from code/)'),
                           ('data', 'files and folders')])
    return _png(fig), W, H


def fig_run():
    """Figure 2: one model run -- the stages of the driver, the modules of
    each, what each leaves."""
    W, H = 17.0, 20.2
    fig, ax = _canvas(W, H)
    xs, ws, xo, wo = 0.35, 10.6, 11.75, 4.9
    ax.text(xs + ws / 2, H - 0.4, 'The run driver, tests/run_lamata_mf6.py: '
            'its stages and their modules', fontsize=7.6, fontweight='bold',
            color=INK, ha='center', va='center')
    ax.text(xo + wo / 2, H - 0.4, 'What the stage leaves', fontsize=7.6,
            fontweight='bold', color=INK, ha='center', va='center')
    stages = [
        ('1  Configuration', 1.42, [
            '`marmites_config.py|reads the .toml, applies --set, validates',
            '`mm_paths.py|the dataset and the workspace of this machine'],
         None),
        ('2  Forcing', 2.1, [
            '`marmites_surface.py|runs the surface model when run.surface '
            'is on:',
            '`   MARMITESsurf/startMARMITESsurface.py, PET_P_INTER.py,',
            '`   MARMITESsurf/rain_seasons.py',
            'otherwise the daily series are read from the dataset'],
         ['`<case>_MMsurf/', 'daily rain, throughfall,',
          'potential ET, LAI']),
        ('3  The model on its grid', 2.45, [
            '`ppMF6/marmites_props.py|layers and their properties',
            '`marmites_meshes.py, marmites_mesh.py|the mesh, the inputs '
            'onto it',
            '`marmites_dem.py|the land surface, from the DEM',
            '`marmites_vector.py, ppMF6/marmites_overlay.py|GIS layers onto '
            'the cells',
            '`ppMF6/marmites_channel.py|the stream cells'],
         ['mesh and overlay caches', '(in the MODFLOW 6 workspace)']),
        ('4  The MODFLOW 6 simulation', 2.45, [
            '`ppMF6/marmites_mf6.py|writes it: unsaturated zone, drains, '
            'solver, time steps',
            '`ppMF6/marmites_sfr.py, marmites_lak.py|streams and ponds',
            '`ppMF6/marmites_evt.py|groundwater ET',
            '`ppMF6/marmites_crr.py|the runoff cascade',
            '`marmites_pest.py|the multipliers, in a calibration run'],
         ['`MF6_ws_<grid>/', 'the MODFLOW 6 input files']),
    ]
    y = H - 0.85
    for title, h, lines, out in stages:
        y -= h
        _box(ax, xs, y, ws, h, 'script', title, lines)
        if out:
            _box(ax, xo, y + (h - 1.3) / 2, wo, 1.3, 'data', None, out)
            _arrow(ax, (xs + ws, y + h / 2), (xo, y + h / 2), color='#7A7A7A',
                   lw=0.7)
        y -= 0.34
        _arrow(ax, (xs + ws / 2, y + 0.34), (xs + ws / 2, y))
    # 5: the coupled march, a loop
    h5 = 3.75
    y -= h5
    _box(ax, xs, y, ws, h5, 'frame2', '5  The coupled march, one day at a '
         'time', ['`marmites_coupler.py|drives MODFLOW 6 through its API'])
    bw, bh, by = 2.36, 1.75, y + 0.72
    inner = [
        ('a  soil step', ['`MARMITESsoil_v3', '`marmites_soil_nb',
                          'infiltration, soil ET,', 'percolation, runoff']),
        ('b  cascade', ['`marmites_crr', 'runoff runs on', 'downslope,',
                        're-infiltrates']),
        ('c  hand over', ['percolation,', 'ET demand, stream',
                          'and pond inflows', 'to MODFLOW 6']),
        ('d  MF6 solves', ['the day; ATS', 'retries shorter', 'steps when',
                           'it fails']),
    ]
    for k, (t, ls) in enumerate(inner):
        bx = xs + 0.3 + k * (bw + 0.22)
        _box(ax, bx, by, bw, bh, 'you', t, ls, fs=6.1, pad=0.12)
        if k:
            _arrow(ax, (bx - 0.22, by + bh / 2), (bx, by + bh / 2))
    x_a = xs + 0.3 + bw / 2
    x_d = xs + 0.3 + 3 * (bw + 0.22) + bw / 2
    ax.plot([x_d, x_d, x_a], [by, by - 0.3, by - 0.3], color='#404040',
            lw=0.9, zorder=3)
    _arrow(ax, (x_a, by - 0.3), (x_a, by))
    ax.text((x_a + x_d) / 2, by - 0.47, 'heads and exfiltration, for the next '
            'day', fontsize=6.2, color=INK, ha='center', va='center')
    _box(ax, xo, y + (h5 - 1.6) / 2, wo, 1.6, 'data', None,
         ['`<name>.hds, <name>.cbc', 'heads and cell budgets,',
          'stream and pond budgets,', 'the listing'])
    _arrow(ax, (xs + ws, y + h5 / 2), (xo, y + h5 / 2), color='#7A7A7A',
           lw=0.7)
    tail = [
        ('6  Save', 1.4, [
            'the end state (heads, soil, unsaturated zone), the coupled '
            'results,',
            'the observation exports: one row per day and point'],
         ['`hi_<state>_*   _coupled_lagged.h5', '`obs_exports/obs_heads.csv,',
          '`   obs_sm, obs_sfr, obs_et']),
        ('7  Figures and tables', 1.75, [
            '`ppMF6/marmites_postprocess.py|maps, time series, budgets',
            '`ppMF6/marmites_cdl_figs.py, marmites_cdl_out.py|inputs, ponds, '
            'run time',
            '`MARMITESplot_v3.py, tests/plot_water_budget.py|older maps, '
            'budget comparison'],
         ['`out_<stamp>_<tag>/_input', '`out_<stamp>_<tag>/_output',
          'PNG figures and CSV tables']),
        ('8  The run\'s record', 1.1, [
            '`marmites_runinfo.py|the configuration as run, the end state'],
         ['`runs/<id>/resolved_config.toml', '`runs/<id>/run_info.json']),
    ]
    for title, h, lines, out in tail:
        _arrow(ax, (xs + ws / 2, y), (xs + ws / 2, y - 0.34))
        y -= 0.34 + h
        _box(ax, xs, y, ws, h, 'script', title, lines)
        oh = max(1.0, min(h, 1.3))
        _box(ax, xo, y + (h - oh) / 2, wo, oh, 'data', None, out)
        _arrow(ax, (xs + ws, y + h / 2), (xo, y + h / 2), color='#7A7A7A',
               lw=0.7)
    return _png(fig), W, H


def fig_pest():
    """Figure 3: the PEST++ suite -- define, set up, launch, follow,
    export, run."""
    W, H = 17.0, 22.3
    fig, ax = _canvas(W, H)
    L, R = 0.3, 16.7
    FW = R - L
    # ---- what it starts from
    y = H - 0.3
    ax.text(L, y - 0.12, 'WHAT A CALIBRATION STARTS FROM', fontsize=6.6,
            fontweight='bold', color=INK, va='center')
    y -= 0.38
    bw = (FW - 0.5) / 3
    hA = 2.1
    y -= hA
    starts = [
        ('A finished run  (page 7)', ['its configuration and its end',
                                      'state: the base of every',
                                      'realisation', '`marmites_runinfo.py']),
        ('The calibration\'s choices  (page 8)', [
            'tabs 1 to 4: base, parameters,', 'observations and weights,',
            'ensemble  →  the [pest] section of', '`configs/<model>.toml']),
        ('The observations', ['the dataset\'s records (page 5);',
                              'MODIS tiles in the RS folder  →',
                              'the catchment\'s ET per composite',
                              '`tools/pest_modis_et.py']),
    ]
    for k, (t, ls) in enumerate(starts):
        _box(ax, L + k * (bw + 0.25), y, bw, hA, 'you' if k < 2 else 'data',
             t, ls)
        _arrow(ax, (L + k * (bw + 0.25) + bw / 2, y), (L + k * (bw + 0.25)
                                                       + bw / 2, y - 0.42))
    # ---- set-up
    y -= 0.42
    hS = 4.85
    y -= hS
    _box(ax, L, y, FW, hS, 'frame', '1  SET-UP    tools/pest_lamata_setup.py'
         '    (page 8, tab 5; about an hour)', [], fs=7.0)
    sw, sh, sy = (FW - 0.9) / 3, 3.85, y + 0.22
    steps = [
        ('org', ['the model as every agent will', 'run it, in',
                 '`<pest_ws>/org/', 'configuration copy, dataset copy,',
                 'mesh and start state, the',
                 'parameter files at their start:', '`mm_mult/   pars/',
                 'and a snapshot of the code in', '`<pest_ws>/code/']),
        ('base', ['one run in org: the test of', 'the whole chain',
                  '`pest_lamata_forward.py', 'then the targets and their',
                  'weights, each group at its', 'share of the starting phi',
                  '`pest_lamata_obs.py', '`org/pest_targets.csv']),
        ('pst', ['pyemu writes', '`<pest_ws>/template/',
                 'pilot points (K, Sy, Ss per', 'layer), one-value parameters,',
                 'template and instruction', 'files, the prior ensemble,',
                 'the control file', '`lamata.pst   forward_run.py']),
    ]
    for k, (t, ls) in enumerate(steps):
        bx = L + 0.25 + k * (sw + 0.2)
        _box(ax, bx, sy, sw, sh, 'pest', t, ls)
        if k:
            _arrow(ax, (bx - 0.2, sy + sh / 2), (bx, sy + sh / 2))
    # ---- launch
    _arrow(ax, (W / 2, y), (W / 2, y - 0.42))
    y -= 0.42
    hL = 6.95
    y -= hL
    _box(ax, L, y, FW, hL, 'frame', '2  LAUNCH    tools/pest_lamata_run.py'
         '    (page 8, tab 5; days)', [], fs=7.0)
    mw, mh = 6.6, 2.8
    my = y + hL - 0.55 - mh
    _box(ax, L + 0.25, my, mw, mh, 'pest', 'The master', [
        '`pestpp-ies   in   <pest_ws>/master/',
        'runs the prior ensemble, then each iteration:',
        'a test of the lambdas on a few realisations,',
        'every realisation with the best one, the update;',
        'drops the realisations with a bad phi, kills',
        'the overdue runs'])
    ax_ = L + 0.25 + mw + 1.5
    aw = R - 0.25 - ax_
    _box(ax, ax_, my, aw, mh, 'pest', 'The agents', [
        '`<pest_ws>/workers/worker_0 ... worker_N',
        'each a copy of template/, with its own dataset',
        'copy and model workspace; up to 20;',
        'more can join a running master:',
        '`pest_lamata_run.py --join N'])
    _arrow(ax, (L + 0.25 + mw, my + mh * 0.66), (ax_, my + mh * 0.66))
    _arrow(ax, (ax_, my + mh * 0.34), (L + 0.25 + mw, my + mh * 0.34))
    ax.text(L + 0.25 + mw + 0.75, my + mh * 0.66 + 0.17, 'parameters',
            fontsize=6.0, color=INK, ha='center')
    ax.text(L + 0.25 + mw + 0.75, my + mh * 0.34 - 0.3, 'model values',
            fontsize=6.0, color=INK, ha='center')
    fh_ = 2.9
    fy = y + 0.22
    _box(ax, L + 0.25, fy, FW - 0.5, fh_, 'script',
         'One forward run, in an agent', [], fs=6.6)
    _arrow(ax, (ax_ + aw / 2, my), (ax_ + aw / 2, fy + fh_))
    fw4 = (FW - 0.5 - 0.3 - 3 * 0.22) / 4
    fsteps = [
        ('a  parameters', ['`forward_run.py', 'pyemu writes the',
                           'model\'s files:', '`mm_mult/*.txt',
                           '`pars/*.csv']),
        ('b  configuration', ['`pest_lamata_forward', 'the one-value',
                              'parameters into the', 'agent\'s .toml']),
        ('c  the model run', ['`run_lamata_mf6.py', 'multipliers applied',
                              'on the run\'s grid:', '`marmites_pest.py']),
        ('d  model values', ['`pest_lamata_obs', 'the observation',
                             'exports  →', '`pest_obs.csv',
                             'results removed']),
    ]
    for k, (t, ls) in enumerate(fsteps):
        bx = L + 0.4 + k * (fw4 + 0.22)
        _box(ax, bx, fy + 0.13, fw4, 2.25, 'data', t, ls, fs=6.1, pad=0.12)
        if k:
            _arrow(ax, (bx - 0.22, fy + 0.13 + 1.12), (bx, fy + 0.13 + 1.12))
    # ---- follow
    _arrow(ax, (W / 2, y), (W / 2, y - 0.42))
    y -= 0.42
    hF = 2.8
    y -= hF
    _box(ax, L, y, FW, hF, 'frame', '3  FOLLOW', [], fs=7.0)
    hw = (FW - 0.7) / 2
    _box(ax, L + 0.25, y + 0.2, hw, 2.0, 'you', 'Page 9, PEST++ progress', [
        '`app/lib/pest_progress.py', 'reads the master\'s record files: where '
        'the run is,', 'the runs sent, back, failed or overdue, phi per '
        'iteration'], fs=6.3)
    _box(ax, L + 0.45 + hw, y + 0.2, hw, 2.0, 'script',
         'A report per iteration', [
             '`tools/pest_lamata_report.py --iter k',
             '`master/report_<k>/', 'phi per group, the ensemble against '
             'the observations,', 'the parameters on their bounds'], fs=6.3)
    # ---- export, run
    _arrow(ax, (W / 2, y), (W / 2, y - 0.42))
    y -= 0.42
    hE = 2.15
    y -= hE
    ew = 9.3
    _box(ax, L, y, ew, hE, 'pest', '4  EXPORT    tools/pest_lamata_export.py'
         '    (page 8, tab 5)', [
             'one realisation of an iteration (base, best or by name)  '
             '→',
             '`<pest_ws>/results/<calibration>_it<N>_<real>/',
             '`configs/<model>_PEST.toml', 'Nothing is run.'], fs=6.6)
    _box(ax, L + ew + 0.55, y, FW - ew - 0.55, hE, 'you',
         '5  RUN THE CALIBRATED MODEL', ['page 7, with <model>_PEST.toml,',
                                         'then page 10 for its figures;',
                                         'its end state: <state>_PEST'],
         fs=6.6)
    _arrow(ax, (L + ew, y + hE / 2), (L + ew + 0.55, y + hE / 2))
    # ---- the jobs
    ax.text(L, y - 0.36, 'Page 8 starts the set-up, the launch and the export '
            'as detached jobs, one at a time:  app/lib/pest_jobs.py  →  '
            'tools/pest_job.py  →', fontsize=6.3, color=INK, va='center')
    ax.text(L, y - 0.7, '<pest_ws>/jobs/<stamp>_<kind>/   (job.log, '
            'status.json, result.json: the end of a job and its exit code)',
            fontsize=6.3, color=INK, va='center')
    return _png(fig), W, H


# ===================================================================== #
#  the document
# ===================================================================== #

def build(figures_dir=None):
    d = Doc()
    f1, f2, f3 = fig_workflow(), fig_run(), fig_pest()
    if figures_dir:
        os.makedirs(figures_dir, exist_ok=True)
        for n, f in (('fig1_workflow', f1), ('fig2_run', f2),
                     ('fig3_pest', f3)):
            with open(os.path.join(figures_dir, n + '.png'), 'wb') as fh:
                fh.write(f[0])

    d.p('MARMITES–MODFLOW 6', 'Title')
    d.p('Guide to the scripts and the workflow', 'Title')
    d.p('What each script does, how they link and where to start, with the '
        'PEST++ calibration suite in detail', 'Subtitle')
    d.p('Branch `MM-MF6_SFR_LAK_CRR`  ·  %s  ·  generated by '
        '`docs/make_scripts_guide.py` (edit the script, not this file)'
        % DATE, 'Subtitle')
    d.p('')
    d.link_list([
        ('s_start', '1  Start here'),
        ('s_where', '2  Where things live'),
        ('s_scripts', '3  The scripts, by role'),
        ('s_run', '4  Inside one model run'),
        ('s_pest', '5  The PEST++ calibration suite'),
        ('s_after', '6  After the calibration: checks against satellites'),
        ('s_cmd', '7  Commands, in one place'),
        ('s_more', '8  Where to read more'),
    ])

    # ------------------------------------------------------------- 1
    d.h(1, '1  Start here', 's_start')
    d.p('MARMITES computes the water balance of the soil, cell by cell and '
        'day by day, and MODFLOW 6 the unsaturated zone, the groundwater, '
        'the streams and the ponds; the two run coupled. You drive '
        'everything from the **front end**, an app in the browser: it edits '
        'one configuration file per model, launches the runs and the '
        'calibration, and shows the results. Every script behind it can '
        'also be started from a command line (section 7).')
    d.p('The order of the work is the order of the app\'s pages:')
    d.steps([
        '**Start the app**: `python code/tools/launch_app.py`, in the '
        'model\'s Python environment. It opens in the browser.',
        '**Home**: say once where things are on this machine (the data, '
        'the run workspace, MODFLOW, the PEST workspace).',
        '**1 Grid**: name the GIS layers (catchment, streams, ponds, DEM), '
        'convert them into the dataset, choose and build the mesh.',
        '**2 to 6**: the surface and the driving forces, the soil, the '
        'unsaturated zone and the groundwater, the state variables (the '
        'observations), the plots. **Validate & save** in the sidebar writes '
        'the configuration, `code/configs/<model>.toml`.',
        '**7 Run**: the validation tab says what is wrong or missing; the '
        'Run tab launches the model and shows its log.',
        '**10 Results**: the figures of each run, maps of any single day, '
        'the comparisons with satellite products.',
        '**8 PEST++ set-up**: define the calibration, set it up, launch it. '
        '**9 PEST++ progress**: follow it.',
        '**8 PEST++ set-up, Export**: write the calibrated model as '
        '`<model>_PEST.toml`, then run it on page 7.',
    ])
    d.figure(f1[0], FIG_W, f1[2] * FIG_W / f1[1],
             'Figure 1. From the data to a calibrated model: each step on an '
             'app page, the script behind it, and what it reads or writes.')
    d.note('**About the names.** The run driver is '
           '`code/tests/run_lamata_mf6.py` and the calibration tools are '
           '`code/tools/pest_lamata_*.py`: they carry La Mata\'s name because '
           'they were built on that case. The driver runs any case (it reads '
           '`paths.case` from the configuration). The calibration tools read '
           'their choices from the configuration too, but still hold a few La '
           'Mata file names (`pest_lamata.toml`, `lamata.pst`): check them '
           'before calibrating a second case.')

    # ------------------------------------------------------------- 2
    d.h(1, '2  Where things live', 's_where')
    d.p('Four places, three of them set on the Home page. The repository '
        'never receives a run\'s output, and a run never opens a shapefile.')
    d.table(['Place', 'Set on', 'What it holds'], [
        ['The repository', '–', [
            '`code/`: every script (section 3)',
            '`code/configs/`: one `.toml` per model, and `paths.local.toml` '
            '(this machine\'s folders, not tracked)',
            '`example/<case>/`: the dataset, the inputs the model reads '
            'directly',
            '`docs/`: the documents']],
        ['Data root', 'Home', [
            '`GIS/`: shapefiles and the DEM, read by the converter only',
            '`RS/`: satellite products (MODIS ET, NDVI scenes, soil '
            'moisture)']],
        ['Run workspace (`WS_ROOT`)', 'Home', [
            '`MF6_ws_<grid>/`: the MODFLOW 6 model, its outputs, the saved '
            'states, the observation exports',
            '`<case>_MMsurf/`: the daily forcing, when the run computes it',
            '`out_<stamp>_<tag>/`: one run\'s figures and tables '
            '(`_input`, `_output`)',
            '`runs/<run id>/`: one run\'s log, status and record']],
        ['PEST workspace', 'Home', [
            '`org/`, `template/`, `master/`, `workers/`, `jobs/`, '
            '`results/`, `code/` (section 5.6). Tens of GB; never inside the '
            'run workspace or the repository']],
    ], [2.4, 1.0, 6.6])

    # ------------------------------------------------------------- 3
    d.h(1, '3  The scripts, by role', 's_scripts')
    d.p('Paths are from `code/`. The tables go from what you start, to the '
        'app, to the model, to what comes after a run.')

    d.h(2, '3.1  What you start')
    d.table(['Script', 'What it does', 'Started from'], [
        ['`tools/launch_app.py`',
         'Starts the front end. On a network drive it first copies the code '
         'to a local disk (a mirror, refreshed at each start) and runs the '
         'app from there; the configurations and the dataset stay in the '
         'repository.', 'command line'],
        ['`tests/run_lamata_mf6.py`',
         '**The run driver.** Builds the MODFLOW 6 model from the '
         'configuration and the dataset, runs the coupled model, saves the '
         'state, draws the figures (section 4). Options: `--config`, '
         '`--set section.key=value` (a one-off override), `--run-tag`, '
         '`--probe`.',
         'page 7, a PEST agent, command line'],
        ['`tools/gis_to_dataset.py`',
         'The converter: the GIS layers the configuration names (catchment, '
         'streams, ponds, soil zones, vegetation, irrigation fields, '
         'observation points, DEM) into the dataset\'s tables. `--dry-run` '
         'previews.', 'page 1, command line'],
        ['`tools/pest_lamata_setup.py`, `pest_lamata_run.py`, '
         '`pest_lamata_report.py`, `pest_lamata_export.py`',
         'The calibration: set-up, launch, report, export (section 5).',
         'page 8, command line'],
        ['`tools/pest_modis_et.py`',
         'MODIS MOD16A2GF tiles into the catchment\'s ET per 8-day '
         'composite; with `--obs` also into the dataset, as the catchment\'s '
         'observed ET.', 'the set-up, command line'],
        ['`tools/validate_grid.py`',
         'Does a run on a mesh mean the same as the run on the structured '
         'grid? Compares their water balances.', 'command line'],
        ['`tools/flags_to_toml.py`',
         'Translates an old command line of the driver into a configuration '
         'file.', 'command line (legacy)'],
    ], [3.0, 5.6, 1.9])

    d.h(2, '3.2  The front end, `app/`')
    d.p('A viewer, an editor and a launcher: it may import the model, the '
        'model never imports it. Runs and calibration jobs are separate '
        'processes, so closing the browser stops nothing.')
    d.table(['File', 'What it does'], [
        ['`Home.py`', 'The first page: the case, and where things are on '
         'this machine (written to `configs/paths.local.toml`).'],
        ['`pages/1_1_-_Grid.py`', 'The GIS layers and their conversion; the '
         'grid (structured, Voronoi mesh or quadtree) and the DEM.'],
        ['`pages/2_2_-_Surface_and_driving_forces.py`', 'Meteo stations, '
         'vegetation types, crops, surface soils; the daily forcing.'],
        ['`pages/3_3_-_Soil.py`', 'Soil zones, thickness, horizons; the '
         'vegetation cover of each cell.'],
        ['`pages/4_4_-_Unsaturated_zone_and_groundwater.py`', 'Layers and '
         'their properties, the unsaturated zone, drains, groundwater ET, '
         'streams, ponds, the runoff cascade, the start state.'],
        ['`pages/5_5_-_State_variables.py`', 'What the run is compared with: '
         'heads, soil moisture, streamflow, actual ET, and the points.'],
        ['`pages/6_6_-_Plots.py`', 'Which figures a run draws.'],
        ['`pages/7_7_-_Run.py`', 'Validation of the configuration, then the '
         'launch and the live log.'],
        ['`pages/8_8_-_PEST_setup.py`', 'The calibration\'s definition (the '
         '`[pest]` section) and its jobs: set up, launch, export.'],
        ['`pages/9_9_-_PEST_progress.py`', 'A calibration\'s progress. '
         'Read-only.'],
        ['`pages/10_10_-_Results.py`', 'The figures of a run by tab; maps of '
         'one day; wetness against a satellite NDVI; satellite soil '
         'moisture.'],
        ['`lib/schema.py`', 'Every field of the panels: label, units, help.'],
        ['`lib/panelui.py`', 'The widgets the pages share.'],
        ['`lib/editor.py`', 'Applies the edits to the configuration, '
         'validates, saves.'],
        ['`lib/checks.py`', 'Everything that can be said about a '
         'configuration before it runs; the Run page lists it and refuses '
         'on it.'],
        ['`lib/runs.py`', 'The run registry: launches a run as a detached '
         'process, its status, its log.'],
        ['`lib/dataset_state.py`', 'Is the dataset still what the GIS '
         'layers say? Starts the converter.'],
        ['`lib/loaders.py`, `lib/results.py`', 'Readers for the pages; which '
         'tab and title each figure gets.'],
        ['`lib/pest_jobs.py`, `lib/pest_progress.py`', 'The calibration\'s '
         'jobs and its progress (section 5).'],
    ], [3.6, 6.4])

    d.h(2, '3.3  Configuration and paths')
    d.table(['Script', 'What it does'], [
        ['`marmites_config.py`', 'The configuration\'s schema: every section '
         'and field of the `.toml`, their defaults and their validation. The '
         'app and the driver both read it.'],
        ['`mm_paths.py`', 'The only machine-specific paths: an `MM_*` '
         'environment variable, else `configs/paths.local.toml` (the Home '
         'page), else a default.'],
        ['`marmites_runinfo.py`', 'What a finished run leaves to build on: '
         'the configuration as it ran and the state it saved. A calibration '
         'starts from both.'],
        ['`mmsurf_config.py`', 'Reads the surface model\'s parameter file as '
         'named records (used by the app and the tests).'],
        ['`marmites_indices.py`', 'The names of the columns of the soil '
         'model\'s result arrays.'],
    ], [3.0, 7.0])

    d.h(2, '3.4  The grid, and the inputs onto it')
    d.table(['Script', 'What it does'], [
        ['`marmites_meshes.py`', 'The mesh producers, chosen by `[grid] '
         'kind`: the structured grid, a quadtree, a Voronoi mesh refined '
         'along the streams and at the ponds.'],
        ['`marmites_gridgen.py`', 'The quadtree, through USGS GRIDGEN.'],
        ['`marmites_mesh.py`', 'Projects the model\'s inputs onto a mesh.'],
        ['`marmites_grid.py`', 'Cell geometry for any grid: area, width.'],
        ['`marmites_dem.py`', 'The land surface from the DEM at its own '
         'resolution, onto the model cells.'],
        ['`marmites_raster.py`, `marmites_vector.py`', 'Rasters sampled at '
         'the cells; points, lines and polygons wrapped onto the cells.'],
        ['`ppMF6/marmites_overlay.py`', 'Polygon layers into cell values: '
         'the majority class, the area mean, the area share (soil zones, '
         'vegetation cover).'],
        ['`ppMF6/marmites_topology.py`', 'Which cells are neighbours on a '
         'mesh.'],
        ['`ppMF6/marmites_channel.py`', 'The stream cells, from the mapped '
         'stream lines.'],
        ['`ppMF6/marmites_props.py`', 'The layers and their properties (K, '
         'Sy, Ss, thickness) from the configuration.'],
    ], [3.4, 6.6])

    d.h(2, '3.5  The model')
    d.table(['Script', 'What it does'], [
        ['`marmites_surface.py`', 'Runs the surface model from the '
         'configuration (when `run.surface` is on) and checks the forcing '
         'files.'],
        ['`MARMITESsurf/startMARMITESsurface.py`, `PET_P_INTER.py`',
         'The surface model: the hourly meteo record into daily rainfall, '
         'interception and throughfall, potential transpiration per '
         'vegetation type, potential evaporation per soil, open-water '
         'evaporation, LAI.'],
        ['`MARMITESsurf/rain_seasons.py`', 'Vegetation seasons started by '
         'the rain: green-up after the first autumn rain, wilting after a '
         'dry spell.'],
        ['`MARMITESsoil/MARMITESsoil_v3.py`', 'The soil model: for a cell '
         'and a day, infiltration, soil evaporation and transpiration, '
         'percolation, runoff, and the groundwater\'s share of ET.'],
        ['`MARMITESsoil/marmites_soil_nb.py`', 'The same soil step, compiled '
         '(numba): what a coupled run uses. A change to the soil step must '
         'be made in both files; a test holds them equal.'],
        ['`ppMF6/marmites_mf6.py`', 'Builds the MODFLOW 6 simulation: grid, '
         'layers, unsaturated zone (UZF), drains, boundaries, solver, '
         'adaptive time steps.'],
        ['`ppMF6/marmites_sfr.py`, `marmites_lak.py`', 'The streams (SFR) '
         'and the ponds (LAK) and their connection.'],
        ['`ppMF6/marmites_evt.py`', 'Groundwater evaporation and '
         'transpiration, as MODFLOW 6 EVT curves.'],
        ['`ppMF6/marmites_crr.py`', 'The runoff cascade: runoff routed to '
         'the downslope cells, where it can re-infiltrate.'],
        ['`marmites_coupler.py`', 'The coupling: for each day, the soil '
         'step, then MODFLOW 6 solved through its API, then the heads and '
         'the exfiltration back to the soil.'],
        ['`marmites_pest.py`', 'The model\'s side of a calibration run: '
         'applies the parameter multipliers on the run\'s own grid.'],
    ], [3.6, 6.4])

    d.h(2, '3.6  After a run')
    d.table(['Script', 'What it does'], [
        ['`ppMF6/marmites_postprocess.py`', 'The figures and tables of a '
         'run: maps, time series at the observation points, water budgets, '
         'Sankey diagrams; and the observation exports the calibration '
         'reads.'],
        ['`ppMF6/marmites_cdl_figs.py`, `marmites_cdl_out.py`', 'Input '
         'figures (parameters, layers, sections through the ponds); pond '
         'recharge, steady-state heads, run time per day.'],
        ['`ppMF6/marmites_rasterise.py`, `marmites_plot.py`', 'A mesh result '
         'as a raster for the map figures; maps for any grid.'],
        ['`MARMITESutilities/MARMITESplot/MARMITESplot_v3.py`', 'The older '
         'figure suite, still drawing the maps and the point time series.'],
        ['`tests/plot_water_budget.py`', 'Water-budget figures, and the '
         'comparison with the older MODFLOW-NWT reference run.'],
        ['`ppMF6/marmites_daymap.py`', 'Maps of one day of a finished run, '
         'for any flux, and the wetness W (Results page).'],
        ['`ppMF6/marmites_rs.py`, `marmites_wetimage.py`', 'Finds and '
         'fetches a Landsat or Sentinel-2 scene; the model\'s wetness '
         'against its NDVI.'],
        ['`ppMF6/marmites_smrs.py`', 'Satellite soil moisture (Copernicus '
         'SWI, SSM) against the model\'s soil saturation.'],
    ], [3.8, 6.2])

    d.h(2, '3.7  Not on the run path')
    d.table(['Scripts', 'What they are'], [
        ['`startMARMITES_v3.py`, `ppMF_FloPy/`', 'The former driver, for '
         'MODFLOW-NWT. Kept as the reference; `ppMF_FloPy` still supplies '
         'the object that holds the model description.'],
        ['`MARMITESutilities/` (the other scripts)', 'Small older helpers '
         '(conversions, plots).'],
        ['`SFR_LAK_CRR/`', 'The Casa de Lobos MODFLOW 6 scripts, kept as the '
         'reference for what was merged into MARMITES.'],
        ['`tests/test_*.py`', 'The test suite (about 1570 tests, half an '
         'hour). `tests/lamata_model.py` builds La Mata for the tests.'],
        ['`tests/diag_*.py`, `tests/diagnose_coupling.py`', 'Diagnostics of '
         'past problems (small MODFLOW 6 models).'],
        ['`docs/make_cookbook.py`, `docs/make_scripts_guide.py`', 'Generate '
         'the cookbook and this guide.'],
    ], [3.6, 6.4])

    # ------------------------------------------------------------- 4
    d.h(1, '4  Inside one model run', 's_run')
    d.p('Page 7 starts the run driver as a separate process (`app/lib/'
        'runs.py`), with the saved configuration. The driver goes through '
        'eight stages (Figure 2). On La Mata (4524 cells, 1949 days) a run '
        'takes about an hour.')
    d.figure(f2[0], FIG_W, f2[2] * FIG_W / f2[1],
             'Figure 2. One model run: the stages of the driver, the modules '
             'each uses, and what each leaves in the run workspace.')
    d.p('Three things are worth knowing:')
    d.bullets([
        '**The outputs in `MF6_ws_<grid>/` are overwritten by the next run** '
        'there. The figures and tables of each run are kept in its own '
        '`out_<stamp>_<tag>/` folder, and its log and record in '
        '`runs/<run id>/`.',
        '**The saved state has a name**, and a run asked to save under the '
        'same name replaces it. A calibration starts from a saved state, so '
        'the PEST tools refuse a state that changed after the run they were '
        'told to start from.',
        '**The forcing is recomputed only when `run.surface` is on.** Off, '
        'the daily series are read from the dataset, whatever the surface '
        'page says.',
    ])

    # ------------------------------------------------------------- 5
    d.h(1, '5  The PEST++ calibration suite', 's_pest')
    d.h(2, '5.1  What it does')
    d.p('The calibration uses **PEST++ IES** (`pestpp-ies`), an iterative '
        'ensemble smoother. It runs the model with many parameter sets at '
        'once (an ensemble of realisations: 100 on La Mata), compares each '
        'run with the observations through one number, the objective '
        'function **phi**, and corrects all the sets together, iteration '
        'after iteration. **pyemu** (Python) writes the PEST files. The '
        'runs are shared among **agents**, each in its own folder, by a '
        '**master**. On La Mata the last calibration had 313 parameters, '
        '8955 observations (7987 fitted), 100 realisations, 3 iterations and '
        '20 agents, and took about two days.')
    d.p('Figure 3 shows the whole suite. The sections after it take the '
        'steps one by one.')
    d.figure(f3[0], FIG_W, f3[2] * FIG_W / f3[1],
             'Figure 3. The PEST++ suite: what a calibration starts from, '
             'the set-up, the launch with its master and agents, one forward '
             'run, the follow-up, the export and the calibrated run.')

    d.h(2, '5.2  Before you start')
    d.bullets([
        '**Home**: the PEST workspace (its own folder, with tens of GB '
        'free), the MODFLOW folder (it must hold `pestpp-ies`), and the '
        'Python executable the agents run the model with.',
        '**A finished run to start from** (page 7): the calibration takes '
        'its configuration as the model to calibrate and its end state as '
        'the start of every realisation.',
        '**The observations** in the dataset (page 5). For the satellite '
        'ET, the MODIS MOD16A2GF tiles in the RS folder, in '
        '`MODIS4<case>/`.',
    ])

    d.h(2, '5.3  The scripts of the suite')
    d.table(['Script', 'Role', 'Reads', 'Writes'], [
        ['`app/pages/8_8_-_PEST_setup.py`',
         'The calibration\'s choices, and the buttons that start its jobs.',
         'the configuration', 'its `[pest]` section'],
        ['`app/lib/pest_jobs.py`, `tools/pest_job.py`',
         'Start the set-up, the launch and the export as detached jobs, one '
         'at a time, and record how each ended.',
         '–', '`jobs/<stamp>_<kind>/`: `job.log`, `status.json`, '
         '`result.json`'],
        ['`tools/pest_lamata_setup.py`',
         'The set-up, in three steps: `org`, `base`, `pst`.',
         'the configuration, the base run, the dataset',
         '`org/`, `template/`, `code/`'],
        ['`tools/pest_lamata_obs.py`',
         'The observations: the targets and their weights at set-up; the '
         'model\'s value of each target after a run.',
         'the observed records; a run\'s observation exports',
         '`pest_targets.csv`; `pest_obs.csv`'],
        ['`tools/pest_modis_et.py`',
         'The catchment\'s ET per MODIS composite.',
         'the MODIS tiles', '`org/modis_et_catchment.csv`'],
        ['`tools/pest_lamata_run.py`',
         'The launch: the master and its agents. `--join N` adds agents to '
         'a running master.',
         '`[pest]` (at launch), `template/`',
         '`master/`, `workers/worker_<k>/`'],
        ['`tools/pest_lamata_forward.py`',
         'One forward run, in an agent: parameters in, model run, model '
         'values out.',
         'the parameter files pyemu wrote', '`pest_obs.csv`'],
        ['`marmites_pest.py`',
         'In the model: the multipliers applied on the run\'s grid.',
         '`mm_mult/*.txt`', '–'],
        ['`app/pages/9_9_-_PEST_progress.py`, `app/lib/pest_progress.py`',
         'The progress of a calibration.',
         'the master\'s record files', 'nothing'],
        ['`tools/pest_lamata_report.py`',
         'Figures and a summary of the prior runs or of an iteration.',
         '`master/`', '`master/report_<k>/`'],
        ['`tools/pest_lamata_export.py`',
         'One realisation as a configuration of its own: the calibrated '
         'model.',
         '`master/`, `template/`',
         '`results/<calibration>_it<N>_<real>/`, '
         '`configs/<model>_PEST.toml`'],
    ], [2.9, 3.3, 2.0, 2.3])

    d.h(2, '5.4  The steps')
    d.h(3, 'Step 0. Define the calibration (page 8, tabs 1 to 4)')
    d.bullets([
        '**Base and workspace**: the finished run every realisation starts '
        'from, and the first day that is fitted (the days before it are a '
        'warm-up: reported, not fitted).',
        '**Parameters**: one row per family and layer or soil horizon: on '
        'or off, starting value, bounds; the spacing of the pilot points.',
        '**Observations and weights**: each group\'s share of the starting '
        'phi and its measurement noise.',
        '**Ensemble and PEST++**: the number of realisations, of iterations '
        'and of agents, and where the ensemble starts (`prior`, `warm`, '
        '`restart`).',
    ])
    d.p('**Validate & save** in the sidebar writes them into the `[pest]` '
        'section of the configuration. The tools read that section; nothing '
        'about a calibration is set anywhere else.')

    d.h(3, 'Step 1. Set up (tab 5, or `pest_lamata_setup.py`)')
    d.p('Three steps, run in a row (`--step all`) or one at a time:')
    d.table(['Step', 'What it does', 'Result'], [
        ['`org`', 'Prepares the model as every agent will run it: a copy of '
         'the configuration (forcing not recomputed, no figures, one cycle '
         'from the base state, multipliers on), a copy of the dataset, the '
         'mesh and the start state, the parameter files at their starting '
         'values, the MODIS series. Takes a snapshot of the code, so that a '
         'calibration of days does not change under a later edit.',
         '`<pest_ws>/org/`, `<pest_ws>/code/`'],
        ['`base`', 'Runs the model once in `org`: the test of the whole '
         'chain. Then builds the targets (one row per observation) and '
         'balances the weights, so that each group carries its share of '
         'the starting phi.',
         '`org/pest_targets.csv`'],
        ['`pst`', 'pyemu writes the PEST files: the pilot points for K, Sy '
         'and Ss per layer, the one-value parameters, the template and '
         'instruction files, the prior ensemble, the control file and the '
         'forward-run script.',
         '`<pest_ws>/template/`'],
    ], [0.9, 6.5, 2.6])
    d.p('It takes about an hour, most of it the base run. The launch is '
        'offered only after a set-up that succeeded.')

    d.h(3, 'Step 2. Launch (tab 5, or `pest_lamata_run.py`)')
    d.p('The launch copies the template into one folder per agent, moves a '
        'previous `master/` aside (as `master_<stamp>`) and starts '
        '`pestpp-ies`. It reads `[pest]` at that moment, so the ensemble '
        'settings can still change after the set-up:')
    d.bullets([
        '**`noptmax`**: `-1` runs the prior ensemble only; `0` makes one '
        'run at the starting values, the check that the whole chain works '
        '(it should give the set-up\'s phi); `n` runs n iterations.',
        '**`start`**: `prior` is the ensemble drawn at set-up; `warm` '
        'starts from a parameter ensemble of an earlier calibration, after '
        'the model changed; `restart` continues an interrupted calibration '
        'of the same model from its master folder.',
        '**`agents`**: up to 20. More can be added while it runs: '
        '`pest_lamata_run.py --join N`.',
    ])
    d.p('`pestpp-ies` then runs the prior ensemble and, in each iteration, '
        'tests its lambda values on a few realisations, runs every '
        'realisation with the best one and updates the parameters. It drops '
        'a realisation whose phi is far above the others and kills a run '
        'that takes much longer than the average.')

    d.h(3, 'Step 3. Follow (page 9, and `pest_lamata_report.py`)')
    d.p('Page 9 reads what the master writes: where the calibration is, the '
        'runs sent, returned, failed or overdue, phi per iteration and per '
        'group, and the end of a failed run\'s log. It changes nothing. '
        'After the prior runs and after each iteration, '
        '`pest_lamata_report.py --iter k` draws a report into '
        '`master/report_<k>/`: phi of each group, the ensemble against the '
        'observations through time, and the parameters that sit on their '
        'bounds.')

    d.h(3, 'Step 4. Export, then run (tab 5, or `pest_lamata_export.py`)')
    d.p('Once an iteration has run, the export writes one of its '
        'realisations as the calibrated model. By default it takes the '
        '`base` realisation of the last iteration; `best` or a name are the '
        'alternatives. It writes the multipliers and the one-value '
        'parameters into `results/<calibration>_it<N>_<real>/` and the '
        'configuration `<model>_PEST.toml` beside the input one. Nothing is '
        'run: choose `<model>_PEST.toml` on page 7 and run it. That run '
        'saves its end state as `<state>_PEST`, never over the base state.')

    d.h(2, '5.5  One forward run, inside an agent')
    d.steps([
        '`pestpp-ies` hands the agent a parameter set. `forward_run.py` '
        '(written by the set-up) lets pyemu turn it into the model\'s files: '
        'the multipliers per cell, interpolated from the pilot points '
        '(`mm_mult/k_l<L>.txt`, `sy_l<L>.txt`, `ss_l<L>.txt`), and the '
        'one-value parameters (`pars/<name>.csv`).',
        '`pest_lamata_forward.py` writes the one-value parameters into the '
        'agent\'s copy of the configuration: the soil horizons\' field '
        'capacity and Ks, the streambed and pond-bed conductances, the '
        'unsaturated zone\'s parameters, the cascade\'s.',
        'It runs the driver in the agent\'s own workspace. While the model '
        'is built, `marmites_pest.py` applies the multipliers on the grid '
        'the model runs on.',
        '`pest_lamata_obs.py` reads the run\'s observation exports '
        '(`obs_heads.csv`, `obs_sm.csv`, `obs_sfr.csv`, `obs_et.csv`) and '
        'writes the model\'s value of every target into `pest_obs.csv`.',
        'The run\'s results folder is removed, so that an agent does not '
        'fill the disk.',
    ])
    d.p('A run that fails still writes `pest_obs.csv`, with a large '
        'sentinel value everywhere, so that PEST++ counts a failed '
        'realisation and goes on. A run that stops converging is ended by '
        'the model\'s own guard (`[run] stall_iterations`, `stall_periods`).')

    d.h(2, '5.6  The PEST workspace')
    d.table(['Folder', 'Written by', 'What it holds'], [
        ['`code/`', 'set-up (`org`)', 'The snapshot of the repository\'s '
         'code the agents run.'],
        ['`org/`', 'set-up (`org`, `base`)', 'The model as an agent runs it '
         '(`pest_lamata.toml`, `example/<case>/`, `ws/`, `mm_mult/`, '
         '`pars/`), the base run\'s exports, `pest_targets.csv`, '
         '`setup.json`.'],
        ['`template/`', 'set-up (`pst`)', 'What every agent is a copy of: '
         '`org` plus the PEST files (`lamata.pst`, the `.tpl` and `.ins` '
         'files, the pilot points, the prior ensemble, `forward_run.py`).'],
        ['`master/`', 'launch', 'The run of `pestpp-ies`: its record '
         '(`lamata.rec`, `lamata.rmr`), phi per iteration '
         '(`lamata.phi.actual.csv`, `lamata.phi.group.csv`), the parameter '
         'and observation ensembles of each iteration '
         '(`lamata.<k>.par.jcb`, `lamata.<k>.obs.jcb`), the reports '
         '(`report_<k>/`).'],
        ['`master_<stamp>/`', 'launch', 'An earlier calibration, moved '
         'aside by the next launch.'],
        ['`workers/worker_<k>/`', 'launch', 'One agent each: a copy of the '
         'template, its own dataset copy and model workspace, its '
         '`model.log`.'],
        ['`jobs/<stamp>_<kind>/`', 'page 8', 'One job each (`setup`, '
         '`launch`, `export`): `job.log`, `status.json`, `result.json`.'],
        ['`results/<calibration>_it<N>_<real>/`', 'export', 'An exported '
         'realisation: its multipliers and parameters, `export.json`.'],
    ], [3.0, 1.7, 5.3])

    d.h(2, '5.7  What is calibrated, against what')
    d.p('The parameter families (page 8, Parameters), each on or off, with '
        'its bounds:')
    d.table(['Family', 'What', 'How'], [
        ['`k`, `sy`, `ss`', 'Hydraulic conductivity, specific yield, '
         'specific storage of each layer', 'pilot points: a value every few '
         'hundred metres, interpolated to the cells'],
        ['`k33`, `sr`, `drn`', 'Kh/Kv, specific retention, drain '
         'conductance of each layer', 'one value per layer'],
        ['`sfc`, `ks`', 'Field capacity and saturated conductivity of each '
         'soil horizon', 'one value per horizon'],
        ['`sfr_rhk`, `lak_bedleak`', 'Streambed conductivity, pond-bed '
         'leakance', 'one value'],
        ['`uzf_eps`, `uzf_vks`', 'The unsaturated zone\'s exponent and '
         'vertical conductivity', 'one value'],
        ['`crr_beta`', 'The share of runoff the cascade passes downslope',
         'one value'],
    ], [2.4, 4.6, 3.0])
    d.p('The observation groups, with La Mata\'s shares of the starting '
        'phi:')
    d.table(['Group', 'What is compared', 'Share'], [
        ['`hds_p0`', 'Measured heads at the main piezometer, near-daily',
         '35 %'],
        ['`hds_c`', 'Heads at the other points: soft data, each weighing a '
         'quarter of a measured head', '(soft)'],
        ['`sm_p0`, `sm_sm`', 'Soil moisture at two probes, two horizons '
         'each', '17.5 % each'],
        ['`q_log`', 'Outlet flow, as its logarithm, daily', '10 %'],
        ['`q_dry`', 'The model\'s flow on the gauge\'s dry days, against '
         'zero', '5 %'],
        ['`q_year`', 'Outlet volume per hydrological year', '5 %'],
        ['`et_rs`', 'MODIS ET over the catchment, per 8-day composite',
         '10 %'],
        ['`et_pt`', 'Measured ET at the points (eddy-covariance tower)',
         '0: reported only'],
    ], [2.2, 6.0, 1.8])

    d.h(2, '5.8  Rules and traps')
    d.bullets([
        '**One PEST job at a time.** A set-up rewrites the template every '
        'agent copies; page 8 refuses a second job while one runs.',
        '**A calibration runs the code as it was at set-up** (the snapshot '
        'in `<pest_ws>/code/`). A later change to the code reaches it only '
        'through a new set-up.',
        '**After a change to the model** (code or configuration): a new '
        'base run on page 7, a new set-up, then a launch. An earlier '
        'ensemble can be reused with `start = warm`.',
        '**Check the chain first**: a launch with `noptmax = 0` makes one '
        'run and should give the set-up\'s phi.',
        '**The PEST workspace is its own folder**: never the run workspace, '
        'never the repository. It grows to tens of GB.',
        '**A job has ended when its `result.json` exists**; `status.json` '
        'only says that it started.',
        '**The export refuses** when the base state or the template changed '
        'after the calibration: the exported model would not be the one '
        'that was calibrated.',
        '**At most 20 agents** in all.',
    ])

    # ------------------------------------------------------------- 6
    d.h(1, '6  After the calibration: checks against satellites', 's_after')
    d.p('These are on the Results page and use the outputs of a finished '
        'run. They compare; they are not part of the calibration.')
    d.table(['Check', 'Script', 'What it gives'], [
        ['Maps of one day', '`ppMF6/marmites_daymap.py`', 'Any flux on any '
         'day, and the wetness W (soil saturation and open water). Only '
         'while the run\'s outputs are still in its workspace.'],
        ['Wetness against NDVI', '`ppMF6/marmites_rs.py`, '
         '`marmites_wetimage.py`', 'Finds and fetches a Landsat or '
         'Sentinel-2 scene; a regression of NDVI on W per cell and '
         'vegetation class; a map of the residuals.'],
        ['Satellite soil moisture', '`ppMF6/marmites_smrs.py`', 'The '
         'Copernicus Soil Water Index or surface soil moisture against the '
         'catchment\'s soil saturation: correlation, bias, by T-value. '
         'Reads product files or the Copernicus Browser\'s CSV exports.'],
    ], [2.2, 3.0, 4.8])

    # ------------------------------------------------------------- 7
    d.h(1, '7  Commands, in one place', 's_cmd')
    d.p('From the repository\'s folder, in the model\'s Python environment. '
        'The app does the same through its pages.')
    d.code([
        '# the app',
        'python code/tools/launch_app.py',
        '',
        '# GIS layers -> dataset (--dry-run to preview)',
        'python code/tools/gis_to_dataset.py --case LaMata',
        '',
        '# a model run; --set overrides one key for this run',
        'python code/tests/run_lamata_mf6.py --config code/configs/'
        'lamata.toml',
        'python code/tests/run_lamata_mf6.py --config code/configs/'
        'lamata.toml --set run.nsp=60 --run-tag try',
        '',
        '# calibration: set-up (all steps, or --step org|base|pst)',
        'python code/tools/pest_lamata_setup.py --config code/configs/'
        'lamata.toml',
        '# launch (settings from [pest]); add agents to a running master',
        'python code/tools/pest_lamata_run.py --config code/configs/'
        'lamata.toml',
        'python code/tools/pest_lamata_run.py --join 4',
        '# report of the prior runs (0) or of iteration k',
        'python code/tools/pest_lamata_report.py --iter 0',
        '# the calibrated model: lamata_PEST.toml, then run it',
        'python code/tools/pest_lamata_export.py --real base',
        'python code/tests/run_lamata_mf6.py --config code/configs/'
        'lamata_PEST.toml',
        '',
        '# MODIS ET over the catchment, also as an observation file',
        'python code/tools/pest_modis_et.py --product MODIS4LaMata --config '
        'code/configs/lamata.toml --obs',
        '',
        '# the tests (all, or one file)',
        'python -m pytest code/tests -q',
        '',
        '# this guide',
        'python docs/make_scripts_guide.py',
    ])

    # ------------------------------------------------------------- 8
    d.h(1, '8  Where to read more', 's_more')
    d.table(['Document', 'What it covers'], [
        ['`CLAUDE.md`', 'The working rules of the project, the '
         'environments, the design decisions in force, the traps already '
         'met.'],
        ['`docs/MM-MF6/CdL_kickoff.md`', 'The state at the end of the La '
         'Mata work and the plan for Casa de Lobos.'],
        ['`docs/MM-MF6/WP7_PESTpp_LaMata_draft.md`', 'The design of the '
         'calibration: parameters, observations, weights.'],
        ['`docs/MM-MF6/WP7_PEST_frontend_plan.md`', 'The PEST pages of the '
         'app, phase by phase.'],
        ['`docs/MM-MF6/MARMITES_SFR_LAK_CRR_analysis.md`', 'The running '
         'notes of the analyses: streams, ponds, cascade, convergence.'],
        ['`docs/MARMITES_x_CdL_merge_cookbook.docx`', 'The plan of the '
         'merge of the Casa de Lobos approach into MARMITES, by work '
         'package.'],
        ['`code/app/README.md`', 'The front end: how to start it, its '
         'layout and its two rules.'],
        ['Each script\'s own first lines', 'Every script opens with a text '
         'that says what it does and why; the tables above are their '
         'summary.'],
    ], [4.0, 6.0])

    d.save(OUT)
    return OUT


if __name__ == '__main__':
    figs = None
    if '--figures' in sys.argv:
        figs = sys.argv[sys.argv.index('--figures') + 1]
    print('written: %s' % build(figs))
