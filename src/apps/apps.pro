TEMPLATE = subdirs
SUBDIRS = builder \
    convert \
    grep

# windows with msvc does not have getopt which
# chemkit-convert and chemkit-grep need
win32-msvc* {
    SUBDIRS -= convert grep
}
