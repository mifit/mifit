import Qt 4.7

Item {
    id: stackArea
    property bool minimized: false

    function updateState() {
        if (stack.visible && stack.size != 0)
            state = minimized ? "minimized" : "shown"
        else
            state = ""
    }

    Component.onCompleted: {
        stack.sizeChanged.connect(updateState)
        stack.visibilityChanged.connect(updateState)
    }

    width: list.width + 20
    anchors.left: root.left
    anchors.leftMargin: 5
    anchors.bottom: root.bottom
    anchors.bottomMargin: 5
    height: (list.height + 7) + stackHeader.height
    opacity: 0
    clip:  true

    Text {
        id: minimizedText
        anchors.left: parent.left
        anchors.leftMargin: 5
        anchors.top: parent.top
        opacity: 0
        color: "white"
        text: "..."
        z: 1
    }
    Rectangle {
        id: stackBox
        color: "gray"
        radius:  5
        opacity: 0.5
        anchors.fill: parent
        z: -1
    }
    Item {
        id: stackHeader
        anchors.top: stackBox.top
        anchors.left: stackBox.left
        anchors.right: stackBox.right
        height: 16
        Rectangle {
            anchors.fill: parent
            color: "blue"
            radius:  5
            opacity: 0.5

        }
        Image {
            id: stackPop
            source: "qrc:/images/popStack.png"
            anchors.topMargin: 2
            anchors.top: parent.top
            anchors.right: stackMinimize.left
            anchors.rightMargin: 5
            MouseArea {
                anchors.fill: parent
                onClicked: { stack.pop() }
            }
        }
        Image {
            id: stackMinimize
            source: (stackArea.minimized ? "qrc:/images/unminimizeStack.png"
                                         : "qrc:/images/minimizeStack.png")
            anchors.topMargin: 2
            anchors.top: parent.top
            anchors.right: stackClear.left
            anchors.rightMargin: 5
            MouseArea {
                anchors.fill: parent
                onClicked: {
                    stackArea.minimized = !stackArea.minimized
                    stackArea.state = (stackArea.minimized) ? "minimized" : "shown";
                }
            }
        }
        Image {
            id: stackClear
            source: "qrc:/images/closeStack.png"
            anchors.topMargin: 2
            anchors.top: parent.top
            anchors.right: parent.right
            anchors.rightMargin: 5
            MouseArea {
                anchors.fill: parent
                onClicked: { stackArea.state = "clear" }
            }
        }
        Text {
            anchors.top: parent.top
            anchors.leftMargin: 5
            anchors.left: parent.left
            visible: stackArea.minimized
            color: "white"
            text: "..."
        }

    }
    Column {
        id: list
        anchors.left: parent.left
        anchors.leftMargin: 10
        anchors.top: stackHeader.bottom
        anchors.topMargin: 3
        Repeater {
            model: stack.stringList
            Text {
                color: "white"
                text: modelData
            }
        }
    }

    states: [
        State {
            name: "minimized"
            PropertyChanges { target: stackArea; height: stackHeader.height; opacity: 1 }
            PropertyChanges { target: minimizedText; opacity: 1 }
        },
        State {
            name: "shown"
            PropertyChanges {
                target: stackArea
                height: (list.height + 7) + stackHeader.height
                opacity: 1
            }
            PropertyChanges { target: minimizedText; opacity: 0 }
        },
        State {
            name:  "clear"
            PropertyChanges { target: stackArea; opacity: 0.0 }
        }
    ]

    transitions: [
        Transition {
            from: ""; to: "minimized"
            SequentialAnimation {
                PropertyAnimation  { target: stackArea; property: "height"; duration: 1 }
                ParallelAnimation {
                    PropertyAnimation  { target: stackArea; property: "opacity"; duration: 300 }
                    PropertyAnimation  { target: minimizedText; property: "opacity"; duration: 300 }
                }
            }
        },
        Transition {
            from: "shown"; to: "minimized"
            PropertyAnimation  { target: stackArea; property: "opacity"; duration: 300 }
            PropertyAnimation  { target: stackArea; property: "height"; duration: 300 }
            PropertyAnimation  { target: minimizedText; property: "opacity"; duration: 300 }
        },
        Transition {
            from: "*"; to: "shown"
            PropertyAnimation  { target: stackArea; property: "opacity"; duration: 300 }
            PropertyAnimation  { target: stackArea; property: "height"; duration: 300 }
            PropertyAnimation  { target: minimizedText; property: "opacity"; duration: 300 }
        },
        Transition {
            from: "*"; to: "clear"
            SequentialAnimation {
                PropertyAnimation  { target: stackArea; property: "opacity"; duration: 300 }
                ScriptAction { script: { stack.clear(); stackArea.state = "" } }
            }
        },
        Transition {
            from: "*"; to: ""
            PropertyAnimation  { target: stackArea; property: "opacity"; duration: 300 }
        }
    ]

}
