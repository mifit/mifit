import Qt 4.7

Item {
    id: stackArea
    property bool minimized: false
    visible: stackVisible && stack.size > 0
    width: list.width + 20
    height: (stackArea.minimized ? 0 : (list.height + 7))
            + stackHeader.height
    anchors.left: root.left
    anchors.leftMargin: 5
    anchors.bottom: root.bottom
    anchors.bottomMargin: 5
    Text {
        id: text
        visible: stackArea.minimized
        anchors.left: parent.left
        anchors.leftMargin: 5
        anchors.top: parent.top
        opacity: 1
        color: "white"
        text: "..."
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
                onClicked: { stackArea.minimized = !stackArea.minimized }
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
                onClicked: { stack.clear() }
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
        visible: !stackArea.minimized
        anchors.left: parent.left
        anchors.leftMargin: 10
        anchors.bottom: parent.bottom
        anchors.bottomMargin: 5
        Repeater {
            model: stack.stringList
            Text {
                id: itemText
                color: "white"
                text: modelData
            }
        }
    }
}
