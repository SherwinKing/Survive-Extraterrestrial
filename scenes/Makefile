.PHONY : all

#n.b. the '-y' sets autoexec scripts to 'on' so that driver expressions will work
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	BLENDER = /Applications/Blender.app/Contents/MacOS/blender -y
else ifeq ($(UNAME_S),Linux)
	BLENDER = ../../blender-2.93.4-linux-x64/blender -y
else
	BLENDER = blender
endif

EXPORT_MESHES=export-meshes.py
EXPORT_WALKMESHES=export-walkmeshes.py
EXPORT_SCENE=export-scene.py

DIST=../dist

all : \
	$(DIST)/space-ground-bomb.pnct \
	$(DIST)/space-ground-bomb.w \
	$(DIST)/space-ground-bomb.scene \

$(DIST)/space-ground-bomb.pnct : space-ground-bomb.blend $(EXPORT_MESHES)
	$(BLENDER) --background --python $(EXPORT_MESHES) -- '$<':Main '$@'

$(DIST)/space-ground-bomb.scene : space-ground-bomb.blend $(EXPORT_SCENE)
	$(BLENDER) --background --python $(EXPORT_SCENE) -- '$<':Main '$@'

echo "$(DIST)/space-ground-bomb.w : space-ground-bomb.blend $(EXPORT_WALKMESHES)
	$(BLENDER) --background --python $(EXPORT_WALKMESHES) -- '$<':WalkMeshes '$@'""
